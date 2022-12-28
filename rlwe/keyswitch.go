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
	LevelSP       map[int]int
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
	levelQ := params.QCount() - 1
	ringQP := params.RingQP()

	for i := 0; i < 6; i++ {
		buff.BuffQP[i] = ringQP.NewPolyLvl(levelQ, levelQ)
	}

	buff.BuffLAQP[0] = ringQP.NewPolyLvl(levelQ, levelQ)
	buff.BuffLAQP[1] = ringQP.NewPolyLvl(levelQ, levelQ)
	buff.BuffLAQP[2] = ringQP.NewPolyLvl(levelQ, levelQ)
	buff.BuffLAQP[3] = ringQP.NewPolyLvl(levelQ, levelQ)

	buff.BuffNTT = params.RingQ().NewPoly()
	buff.BuffInvNTT = params.RingQ().NewPoly()

	buff.BuffDecompQP = make([]PolyQP, params.Beta())
	for i := 0; i < params.Beta(); i++ {
		buff.BuffDecompQP[i] = ringQP.NewPolyLvl(levelQ, levelQ)
	}

	return buff
}

// NewKeySwitcher creates a new KeySwitcher.
func NewKeySwitcher(params Parameters) *KeySwitcher {

	ks := new(KeySwitcher)
	ks.Parameters = &params

	// Generate necessary rings for level aware gadget decomposition
	pCount := params.PCount()
	qCount := params.QCount()

	if qCount%pCount != 0 {
		panic("pCount should divide qCount!!")
	}

	ks.RingPk = make([]*ring.Ring, params.Beta())
	ks.RingQPk = make([]RingQP, params.Beta())
	ks.PkDivP = make([]*ring.Poly, params.Beta())

	ringQ := params.RingQ()

	for i := range ks.RingPk {
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

	// Generate BasisExtender and Decomposer for each special modulus
	ks.BasisExtender = make([]*ring.BasisExtender, params.Beta())
	ks.Decomposer = make([]*ring.Decomposer, params.Beta())
	for i := range ks.BasisExtender {
		ks.BasisExtender[i] = ring.NewBasisExtender(ringQ, ks.RingPk[i])
		ks.Decomposer[i] = ring.NewDecomposer(ringQ, ks.RingPk[i])
	}

	ks.keySwitcherBuffer = newKeySwitcherBuffer(params)

	// Generate proper special modulus index for each level
	// TODO: is it optimal?
	ks.LevelSP = make(map[int]int)
	for levelQ := 0; levelQ < qCount; levelQ++ {
		min_cost := 987654321

		for i := 1; i <= params.Beta()/2; i++ {
			levelSP := pCount*i - 1
			if levelQ+levelSP+2 > qCount+pCount {
				break
			}

			decompSize := int(math.Ceil(float64(levelQ+1) / float64(levelSP+1)))
			cost := (decompSize + 2) * (levelQ + levelSP + 2)

			if cost < min_cost {
				min_cost = cost
				ks.LevelSP[levelQ] = levelSP
			}
		}
	}

	return ks
}

// ShallowCopy creates a copy of a KeySwitcher, only reallocating the memory Buff.
func (ks *KeySwitcher) ShallowCopy() *KeySwitcher {
	return NewKeySwitcher(*ks.Parameters)
}

/*
func (ks *KeySwitcher) LevelSP(levelQ int) int {
	return ks.LevelSP[levelQ]
}
*/

// ExtendSpecialModulus rearrages polyQP so that it can have additional special modulus
func (ks *KeySwitcher) ExtendSpecialModulus(levelQ int, polyQPIn, polyQPOut PolyQP) {
	levelSP := ks.LevelSP[levelQ]
	pCount := ks.PCount()
	qCount := ks.QCount()

	for i := 0; i <= levelQ; i++ {
		polyQPOut.Q.Coeffs[i] = polyQPIn.Q.Coeffs[i]
	}

	k := levelSP / pCount
	for i := 0; i < k*pCount; i++ {
		polyQPOut.P.Coeffs[i] = polyQPIn.Q.Coeffs[qCount-k*pCount+i]
	}
	for i := k * pCount; i <= levelSP; i++ {
		polyQPOut.P.Coeffs[i] = polyQPIn.P.Coeffs[i-k*pCount]
	}
}

// SwitchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
// Will return the result in the same NTT domain as the input cx.
func (ks *KeySwitcher) SwitchKeysInPlace(levelQ int, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {
	ks.SwitchKeysInPlaceNoModDown(levelQ, cx, evakey, p0, ks.BuffQP[1].P, p1, ks.BuffQP[2].P)

	levelSP := ks.LevelSP[levelQ]
	k := levelSP / ks.PCount()

	if cx.IsNTT {
		ks.BasisExtender[k].ModDownQPtoQNTT(levelQ, levelSP, p0, ks.BuffQP[1].P, p0)
		ks.BasisExtender[k].ModDownQPtoQNTT(levelQ, levelSP, p1, ks.BuffQP[2].P, p1)
	} else {

		ks.RingQ().InvNTTLazyLvl(levelQ, p0, p0)
		ks.RingQ().InvNTTLazyLvl(levelQ, p1, p1)
		ks.RingPk[k].InvNTTLazyLvl(levelSP, ks.BuffQP[1].P, ks.BuffQP[1].P)
		ks.RingPk[k].InvNTTLazyLvl(levelSP, ks.BuffQP[2].P, ks.BuffQP[2].P)

		ks.BasisExtender[k].ModDownQPtoQ(levelQ, levelSP, p0, ks.BuffQP[1].P, p0)
		ks.BasisExtender[k].ModDownQPtoQ(levelQ, levelSP, p1, ks.BuffQP[2].P, p1)
	}
}

// DecomposeNTT applies the full RNS basis decomposition for all q_alpha_i on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// BuffDecompQ and BuffDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (ks *KeySwitcher) DecomposeNTT(levelQ, levelSP, alpha int, c2 *ring.Poly, BuffDecomp []PolyQP) {

	ringQ := ks.RingQ()
	k := levelSP / ks.PCount()

	polyNTT := ks.BuffNTT
	polyInvNTT := ks.BuffInvNTT

	if c2.IsNTT {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, c2, ks.PkDivP[k], polyNTT)
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, c2, ks.PkDivP[k], polyInvNTT)
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	decompSize := int(math.Ceil(float64(levelQ+1) / float64(levelSP+1)))

	for i := 0; i < decompSize; i++ {
		ks.DecomposeSingleNTT(levelQ, levelSP, alpha, i, polyNTT, polyInvNTT, BuffDecomp[i].Q, BuffDecomp[i].P)
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (ks *KeySwitcher) DecomposeSingleNTT(levelQ, levelSP, alpha, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	k := levelSP / ks.PCount()
	ringQ := ks.RingQ()
	ringPk := ks.RingPk[k]

	ks.Decomposer[k].DecomposeAndSplit(levelQ, levelSP, alpha, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * (levelSP + 1)
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
	ringPk.NTTLazyLvl(levelSP, c2QiP, c2QiP)
}

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// Buff2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// Buff3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (ks *KeySwitcher) SwitchKeysInPlaceNoModDown(levelQ int, cx *ring.Poly, evakey *SwitchingKey, c0Q, c0P, c1Q, c1P *ring.Poly) {

	levelSP := ks.LevelSP[levelQ]
	k := levelSP / ks.PCount()

	ringQ := ks.RingQ()
	ringPk := ks.RingPk[k]
	ringQPk := ks.RingQPk[k]

	c2QP := ks.BuffQP[0]

	cxNTT := ks.BuffNTT
	cxInvNTT := ks.BuffInvNTT
	if cx.IsNTT {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, cx, ks.PkDivP[k], cxNTT)
		ringQ.InvNTTLvl(levelQ, cxNTT, cxInvNTT)
	} else {
		ringQ.MulCoeffsMontgomeryLvl(levelQ, cx, ks.PkDivP[k], cxInvNTT)
		ringQ.NTTLvl(levelQ, cxInvNTT, cxNTT)
	}

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	OverFlowMargin := ks.QiOverflowMargin(ks.QCount()-1) >> 1
	if OverFlowMargin > (ks.PiOverflowMargin(ks.PCount()-1) >> 1) {
		OverFlowMargin = ks.PiOverflowMargin(ks.PCount()-1) >> 1
	}

	decompSize := int(math.Ceil(float64(levelQ+1) / float64(levelSP+1)))

	reduce := 0
	// Key switching with CRT decomposition for the Qi
	for i := 0; i < decompSize; i++ {
		ks.DecomposeSingleNTT(levelQ, levelSP, levelSP+1, i, cxNTT, cxInvNTT, c2QP.Q, c2QP.P)

		ks.ExtendSpecialModulus(levelQ, evakey.Value[levelSP][i][0], ks.BuffLAQP[0])
		ks.ExtendSpecialModulus(levelQ, evakey.Value[levelSP][i][1], ks.BuffLAQP[1])

		if i == 0 {
			ringQPk.MulCoeffsMontgomeryConstantLvl(levelQ, levelSP, ks.BuffLAQP[0], c2QP, c0QP)
			ringQPk.MulCoeffsMontgomeryConstantLvl(levelQ, levelSP, ks.BuffLAQP[1], c2QP, c1QP)
		} else {
			ringQPk.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelSP, ks.BuffLAQP[0], c2QP, c0QP)
			ringQPk.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelSP, ks.BuffLAQP[1], c2QP, c1QP)
		}

		if reduce%OverFlowMargin == OverFlowMargin-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)

			ringPk.ReduceLvl(levelSP, c0QP.P, c0QP.P)
			ringPk.ReduceLvl(levelSP, c1QP.P, c1QP.P)
		}

		reduce++
	}

	if reduce%OverFlowMargin != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)

		ringPk.ReduceLvl(levelSP, c0QP.P, c0QP.P)
		ringPk.ReduceLvl(levelSP, c1QP.P, c1QP.P)
	}

}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod Q
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod Q
func (ks *KeySwitcher) KeyswitchHoisted(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	ks.KeyswitchHoistedNoModDown(levelQ, BuffDecompQP, evakey, c0Q, c1Q, c0P, c1P)

	levelSP := ks.LevelSP[levelQ]
	k := levelSP / ks.PCount()

	// Computes c0Q = c0Q/c0P and c1Q = c1Q/c1P
	ks.BasisExtender[k].ModDownQPtoQNTT(levelQ, levelSP, c0Q, c0P, c0Q)
	ks.BasisExtender[k].ModDownQPtoQNTT(levelQ, levelSP, c1Q, c1P, c1Q)
}

// KeyswitchHoistedNoModDown applies the key-switch to the decomposed polynomial c2 mod QP (BuffDecompQ and BuffDecompP)
//
// Buff2 = dot(BuffDecompQ||BuffDecompP * evakey[0]) mod QP
// Buff3 = dot(BuffDecompQ||BuffDecompP * evakey[1]) mod QP
func (ks *KeySwitcher) KeyswitchHoistedNoModDown(levelQ int, BuffDecompQP []PolyQP, evakey *SwitchingKey, c0Q, c1Q, c0P, c1P *ring.Poly) {

	levelSP := ks.LevelSP[levelQ]
	k := levelSP / ks.PCount()

	ringQ := ks.RingQ()
	ringPk := ks.RingPk[k]
	ringQPk := ks.RingQPk[k]

	c0QP := PolyQP{c0Q, c0P}
	c1QP := PolyQP{c1Q, c1P}

	OverFlowMargin := ks.QiOverflowMargin(ks.QCount()-1) >> 1
	if OverFlowMargin > (ks.PiOverflowMargin(ks.PCount()-1) >> 1) {
		OverFlowMargin = ks.PiOverflowMargin(ks.PCount()-1) >> 1
	}

	decompSize := int(math.Ceil(float64(levelQ+1) / float64(levelSP+1)))

	// Key switching with CRT decomposition for the Qi
	reduce := 0
	for i := 0; i < decompSize; i++ {

		ks.ExtendSpecialModulus(levelQ, evakey.Value[levelSP][i][0], ks.BuffLAQP[0])
		ks.ExtendSpecialModulus(levelQ, evakey.Value[levelSP][i][1], ks.BuffLAQP[1])

		if i == 0 {
			ringQPk.MulCoeffsMontgomeryConstantLvl(levelQ, levelSP, ks.BuffLAQP[0], BuffDecompQP[i], c0QP)
			ringQPk.MulCoeffsMontgomeryConstantLvl(levelQ, levelSP, ks.BuffLAQP[1], BuffDecompQP[i], c1QP)
		} else {
			ringQPk.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelSP, ks.BuffLAQP[0], BuffDecompQP[i], c0QP)
			ringQPk.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, levelSP, ks.BuffLAQP[1], BuffDecompQP[i], c1QP)
		}

		if reduce%OverFlowMargin == OverFlowMargin-1 {
			ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
			ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)

			ringPk.ReduceLvl(levelSP, c0QP.P, c0QP.P)
			ringPk.ReduceLvl(levelSP, c1QP.P, c1QP.P)
		}

		reduce++

	}

	if reduce%OverFlowMargin != 0 {
		ringQ.ReduceLvl(levelQ, c0QP.Q, c0QP.Q)
		ringQ.ReduceLvl(levelQ, c1QP.Q, c1QP.Q)

		ringPk.ReduceLvl(levelSP, c0QP.P, c0QP.P)
		ringPk.ReduceLvl(levelSP, c1QP.P, c1QP.P)
	}
}
