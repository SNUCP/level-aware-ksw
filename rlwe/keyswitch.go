package rlwe

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"math"
)

// KeySwitcher is a struct for RLWE key-switching.
type KeySwitcher struct {
	*Parameters
	*keySwitcherBuffer

	BasisExtender []*ring.BasisExtender
	Decomposer    []*ring.Decomposer
	RingPk        []*ring.Ring
	RingQPk       []RingQP
	SPIndex       []int
	PkDivP        []*ring.Poly
}

type keySwitcherBuffer struct {
	// BuffQ[0]/BuffP[0] : on the fly decomp(c2)
	// BuffQ[1-5]/BuffP[1-5] : available
	BuffQP       [6]PolyQP
	BuffLAQP     [4]PolyQP
	BuffNTT      *ring.Poly
	BuffInvNTT   *ring.Poly
	BuffDecompQP []PolyQP // Memory Buff for the basis extension in hoisting
}

func newKeySwitcherBuffer(params Parameters) *keySwitcherBuffer {

	buff := new(keySwitcherBuffer)
	beta := params.Beta()
	levelQ := params.QCount() - 1
	levelP := params.PCount()*(beta+1) - 1
	ringQP := params.RingQP()

	for i := 0; i < 6; i++ {
		buff.BuffQP[i] = ringQP.NewPolyLvl(levelQ, levelP)
	}

	buff.BuffLAQP[0] = ringQP.NewPolyLvl(levelQ, levelP)
	buff.BuffLAQP[1] = ringQP.NewPolyLvl(levelQ, levelP)
	buff.BuffLAQP[2] = ringQP.NewPolyLvl(levelQ, levelP)
	buff.BuffLAQP[3] = ringQP.NewPolyLvl(levelQ, levelP)

	buff.BuffNTT = params.RingQ().NewPoly()
	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]PolyQP, beta)
	for i := 0; i < beta; i++ {
		buff.BuffDecompQP[i] = ringQP.NewPolyLvl(levelQ, levelP)
	}

	return buff
}

// NewKeySwitcher creates a new KeySwitcher.
func NewKeySwitcher(params Parameters) *KeySwitcher {

	if params.QCount()%params.PCount() != 0 {
		panic("NewKeySwitcher: pcount does not divide qcount, cannot apply level-aware gadget")
	}

	ks := new(KeySwitcher)
	ks.Parameters = &params

	// Generate necessary rings for level aware gadget decomposition
	pCount := params.PCount()
	qCount := params.QCount()
	maxSPIdx := (pCount + qCount) / pCount

	ks.RingPk = make([]*ring.Ring, maxSPIdx)
	ks.RingQPk = make([]RingQP, maxSPIdx)
	ks.PkDivP = make([]*ring.Poly, maxSPIdx)

	ringQ := params.RingQ()

	for i := 0; i < maxSPIdx; i++ {
		pi := append([]uint64{}, params.qi[qCount-i*pCount:]...)
		pi = append(pi, params.pi...)

		ks.RingPk[i], _ = ring.NewRingFromType(1<<params.logN, pi, params.ringType)
		ks.RingQPk[i] = RingQP{RingQ: ringQ, RingP: ks.RingPk[i]}
		ks.PkDivP[i] = ringQ.NewPoly()

		if i == 0 {
			ringQ.AddScalar(ks.PkDivP[i], 1, ks.PkDivP[i])
			ringQ.MForm(ks.PkDivP[i], ks.PkDivP[i])
		} else {
			pkDivP := ks.RingPk[i].ModulusAtLevel[i*pCount-1]
			ringQ.AddScalarBigint(ks.PkDivP[i], pkDivP, ks.PkDivP[i])
			ringQ.MForm(ks.PkDivP[i], ks.PkDivP[i])
		}
	}

	// Generate proper special modulus index for each level
	// TODO: is it optimal?
	ks.SPIndex = make([]int, qCount)
	for level := 0; level < qCount; level++ {
		// ks.SPIndex[level] = 15 // BENCHMARK: Changed Here
		ks.SPIndex[level] = (qCount - 1 - level) / pCount
	}

	// Generate BasisExtender and Decomposer for each special modulus
	ks.BasisExtender = make([]*ring.BasisExtender, maxSPIdx)
	ks.Decomposer = make([]*ring.Decomposer, maxSPIdx)
	for i := 0; i < maxSPIdx; i++ {
		ks.BasisExtender[i] = ring.NewBasisExtender(ringQ, ks.RingPk[i])
		ks.Decomposer[i] = ring.NewDecomposer(ringQ, ks.RingPk[i])
	}

	ks.keySwitcherBuffer = newKeySwitcherBuffer(params)

	return ks
}

// ShallowCopy creates a copy of a KeySwitcher, only reallocating the memory Buff.
func (ks *KeySwitcher) ShallowCopy() *KeySwitcher {
	return NewKeySwitcher(*ks.Parameters)
}

// ExtendSpecialModulus rearrages polyQP so that it can have additional special modulus
func (ks *KeySwitcher) ExtendSpecialModulus(levelQ int, polyQPIn, polyQPOut PolyQP) {
	sp := ks.SPIndex[levelQ]
	pCount := ks.PCount()
	qCount := ks.QCount()

	for l := 0; l < levelQ+1; l++ {
		polyQPOut.Q.Coeffs[l] = polyQPIn.Q.Coeffs[l]
	}

	for l := 0; l < sp*pCount; l++ {
		polyQPOut.P.Coeffs[l] = polyQPIn.Q.Coeffs[qCount-sp*pCount+l]
	}

	for l := sp * pCount; l < (sp+1)*pCount; l++ {
		polyQPOut.P.Coeffs[l] = polyQPIn.P.Coeffs[l-sp*pCount]
	}
}

// LevelPk returns proper special modulus level given levelQ
func (ks *KeySwitcher) LevelPk(levelQ int) int {
	sp := ks.SPIndex[levelQ]
	return (sp+1)*ks.PCount() - 1
}

// SwitchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
// Will return the result in the same NTT domain as the input cx.
func (ks *KeySwitcher) SwitchKeysInPlace(levelQ int, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {
	ks.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, ks.BuffQP[1].P, p1, ks.BuffQP[2].P)

	sp := ks.SPIndex[levelQ]
	levelP := ks.LevelPk(levelQ)

	if cx.IsNTT {
		ks.BasisExtender[sp].ModDownQPtoQNTT(levelQ, levelP, p0, ks.BuffQP[1].P, p0)
		ks.BasisExtender[sp].ModDownQPtoQNTT(levelQ, levelP, p1, ks.BuffQP[2].P, p1)
	} else {

		ks.RingQ().InvNTTLazyLvl(levelQ, p0, p0)
		ks.RingQ().InvNTTLazyLvl(levelQ, p1, p1)
		ks.RingPk[sp].InvNTTLazyLvl(levelP, ks.BuffQP[1].P, ks.BuffQP[1].P)
		ks.RingPk[sp].InvNTTLazyLvl(levelP, ks.BuffQP[2].P, ks.BuffQP[2].P)

		ks.BasisExtender[sp].ModDownQPtoQ(levelQ, levelP, p0, ks.BuffQP[1].P, p0)
		ks.BasisExtender[sp].ModDownQPtoQ(levelQ, levelP, p1, ks.BuffQP[2].P, p1)
	}
}

// DecomposeNTT applies the full RNS basis decomposition for all q_alpha_i on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// BuffDecompQ and BuffDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (ks *KeySwitcher) DecomposeNTT(levelQ, levelP, alpha int, c2 *ring.Poly, BuffDecomp []PolyQP) {

	ringQ := ks.RingQ()
	sp := ks.SPIndex[levelQ]

	polyNTT := ks.BuffNTT
	polyInvNTT := ks.BuffInvNTT

	if c2.IsNTT {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, c2, ks.PkDivP[sp], polyNTT)
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, c2, ks.PkDivP[sp], polyInvNTT)
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	beta := int(math.Ceil(float64(levelQ+1) / float64(ks.PCount())))

	for i := 0; i < beta; i++ {
		if i%(sp+1) == 0 {
			ks.DecomposeSingleNTT(levelQ, levelP, alpha, i/(sp+1), polyNTT, polyInvNTT, BuffDecomp[i/(sp+1)].Q, BuffDecomp[i/(sp+1)].P)
		}
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (ks *KeySwitcher) DecomposeSingleNTT(levelQ, levelP, alpha, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	sp := ks.SPIndex[levelQ]
	ringQ := ks.RingQ()
	ringP := ks.RingPk[sp]

	ks.Decomposer[sp].DecomposeAndSplit(levelQ, levelP, alpha, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * (levelP + 1)
	p0idxed := p0idxst + 1

	// c2_qi = cx mod qi mod qi
	for x := 0; x < levelQ+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ringQ.NTTSingleLazy(x, c2QiQ.Coeffs[x], c2QiQ.Coeffs[x])
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazyLvl(levelP, c2QiP, c2QiP)
}

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// Buff2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// Buff3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (ks *KeySwitcher) SwitchKeysInPlaceNoModDown(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	sp := ks.SPIndex[levelQ]
	ringQ := ks.RingQ()
	ringP := ks.RingPk[sp]
	ringQP := ks.RingQPk[sp]

	c2QP := ks.BuffQP[0]

	cxNTT := ks.BuffNTT
	cxInvNTT := ks.BuffInvNTT
	if cx.IsNTT {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, cx, ks.PkDivP[sp], cxNTT)
		ringQ.InvNTTLvl(levelQ, cxNTT, cxInvNTT)
	} else {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, cx, ks.PkDivP[sp], cxInvNTT)
		ringQ.NTTLvl(levelQ, cxInvNTT, cxNTT)
	}

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	levelP := ks.LevelPk(levelQ)
	beta := int(math.Ceil(float64(levelQ+1) / float64(ks.PCount())))

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(ks.PCount()-1) >> 1

	reduce := 0
	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i += (sp + 1) {

		ks.DecomposeSingleNTT(levelQ, levelP, levelP+1, i/(sp+1), cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		if (sp == 0) || (i+1 == beta) {
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][0], ks.BuffLAQP[0])
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][1], ks.BuffLAQP[1])

			if i == 0 {
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[0], c2QP, c0QP)
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[1], c2QP, c1QP)
			} else {
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], c2QP, c0QP)
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[1], c2QP, c1QP)
			}

		} else {
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][0], ks.BuffLAQP[0])
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i+1][0], ks.BuffLAQP[1])
			ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], ks.BuffLAQP[1], ks.BuffLAQP[2])

			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][1], ks.BuffLAQP[0])
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i+1][1], ks.BuffLAQP[1])
			ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], ks.BuffLAQP[1], ks.BuffLAQP[3])

			for j := 2; j < (sp+1) && (i+j < beta); j++ {
				ks.ExtendSpecialModulus(levelQ, evakey.Value[i+j][0], ks.BuffLAQP[0])
				ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], ks.BuffLAQP[2], ks.BuffLAQP[2])

				ks.ExtendSpecialModulus(levelQ, evakey.Value[i+j][1], ks.BuffLAQP[1])
				ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[1], ks.BuffLAQP[3], ks.BuffLAQP[3])
			}

			if i == 0 {
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[2], c2QP, c0QP)
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[3], c2QP, c1QP)
			} else {
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[2], c2QP, c0QP)
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[3], c2QP, c1QP)
			}
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}

}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod Q
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod Q
func (ks *KeySwitcher) KeyswitchHoisted(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ks.KeyswitchHoistedNoModDown(levelQ, BuffDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	sp := ks.SPIndex[levelQ]
	levelP := ks.LevelPk(levelQ)

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	ks.BasisExtender[sp].ModDownQPtoQNTT(levelQ, levelP, c0Q, c0P, c0Q)
	ks.BasisExtender[sp].ModDownQPtoQNTT(levelQ, levelP, c1Q, c1P, c1Q)
}

// KeyswitchHoistedNoModDown applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod QP
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod QP
func (ks *KeySwitcher) KeyswitchHoistedNoModDown(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	sp := ks.SPIndex[levelQ]
	ringQ := ks.RingQ()
	ringP := ks.RingPk[sp]
	ringQP := ks.RingQPk[sp]

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	levelP := ks.LevelPk(levelQ)
	beta := int(math.Ceil(float64(levelQ+1) / float64(ks.PCount())))

	QiOverF := ks.Parameters.QiOverflowMargin(levelQ) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin(ks.PCount()-1) >> 1

	// Key switching with CRT decomposition for the Qi
	reduce := 0
	for i := 0; i < beta; i += (sp + 1) {

		if (sp == 0) || (i+1 == beta) {

			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][0], ks.BuffLAQP[0])
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][1], ks.BuffLAQP[1])

			if i == 0 {
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[0], BuffDecompQP[i/(sp+1)], c0QP)
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[1], BuffDecompQP[i/(sp+1)], c1QP)
			} else {
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], BuffDecompQP[i/(sp+1)], c0QP)
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[1], BuffDecompQP[i/(sp+1)], c1QP)
			}
		} else {
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][0], ks.BuffLAQP[0])
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i+1][0], ks.BuffLAQP[1])
			ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], ks.BuffLAQP[1], ks.BuffLAQP[2])

			ks.ExtendSpecialModulus(levelQ, evakey.Value[i][1], ks.BuffLAQP[0])
			ks.ExtendSpecialModulus(levelQ, evakey.Value[i+1][1], ks.BuffLAQP[1])
			ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], ks.BuffLAQP[1], ks.BuffLAQP[3])

			for j := 2; j < (sp+1) && (i+j < beta); j++ {
				ks.ExtendSpecialModulus(levelQ, evakey.Value[i+j][0], ks.BuffLAQP[0])
				ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[0], ks.BuffLAQP[2], ks.BuffLAQP[2])

				ks.ExtendSpecialModulus(levelQ, evakey.Value[i+j][1], ks.BuffLAQP[1])
				ringQP.AddNoModLvl(levelQ, levelP, ks.BuffLAQP[1], ks.BuffLAQP[3], ks.BuffLAQP[3])
			}

			if i == 0 {
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[2], BuffDecompQP[i/(sp+1)], c0QP)
				ringQP.MulCoeffsMontgomeryConstantLvl(levelQ, levelP, ks.BuffLAQP[3], BuffDecompQP[i/(sp+1)], c1QP)
			} else {
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[2], BuffDecompQP[i/(sp+1)], c0QP)
				ringQP.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelP, ks.BuffLAQP[3], BuffDecompQP[i/(sp+1)], c1QP)
			}
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
			ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
		}

		reduce++

	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)
	}

	if reduce%PiOverF != 0 {
		ringP.ReduceLvl(levelP, c0QP.P, c0QP.P)
		ringP.ReduceLvl(levelP, c1QP.P, c1QP.P)
	}
}
