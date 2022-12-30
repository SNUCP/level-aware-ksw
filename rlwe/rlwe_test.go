package rlwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ring"

	"github.com/stretchr/testify/require"
)

var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
//var TestParams = []ParametersLiteral{TestPN12QP109, TestPN13QP218, TestPN14QP438, TestPN15QP880, TestPN16QP240, TestPN17QP360}

var rlweTestParam = ParametersLiteral{
	LogN: 16,
	Q: []uint64{
		// 44 x 40
		0xffff480001, 0xffff420001, 0xffff340001, 0xfffeb60001,
		0xfffeb00001, 0xfffe9e0001, 0xfffe860001, 0xfffe680001,
		0xfffe620001, 0xfffe4a0001, 0xfffe2c0001, 0xfffe100001,
		0xfffd800001, 0xfffd720001, 0xfffd6e0001, 0xfffd5a0001,

		0xfffd3e0001, 0xfffd260001, 0xfffd080001, 0xfffcfa0001,
		0xfffcf60001, 0xfffcc60001, 0xfffca00001, 0xfffc940001,
		0xfffc880001, 0xfffc6a0001, 0xfffc640001, 0xfffc600001,
		0xfffc540001, 0xfffc360001, 0xfffc1e0001, 0xfffbf40001,

		0xfffbdc0001, 0xfffbb80001, 0xfffba60001, 0xfffba00001,
		0xfffb5e0001, 0xfffb340001, 0xfffb1a0001, //0xfffb0e0001,
	},
	P: []uint64{
		// 60
		0x1000000000ce0001,
		//0xffff8a0001, 0xffff820001, 0xffff780001, //0xffff580001,
	},
	RingType: ring.Standard,
}

var TestParams = []ParametersLiteral{rlweTestParam}

func testString(params Parameters, opname string) string {
	return fmt.Sprintf("%s/logN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestRLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams[:] {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		kgen := NewKeyGenerator(params)

		for _, testSet := range []func(kgen KeyGenerator, t *testing.T){
			testGenKeyPair,
			testSwitchKeyGen,
			testEncryptor,
			testDecryptor,
			testKeySwitcher,
		} {
			testSet(kgen, t)
			runtime.GC()
		}
	}
}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients
func log2OfInnerSum(level int, ringQ *ring.Ring, poly *ring.Poly) (logSum int) {
	sumRNS := make([]uint64, level+1)
	var sum uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		qiHalf := qi >> 1
		coeffs := poly.Coeffs[i]
		sum = 0

		for j := 0; j < ringQ.N; j++ {

			v := coeffs[j]

			if v >= qiHalf {
				sum = ring.CRed(sum+qi-v, qi)
			} else {
				sum = ring.CRed(sum+v, qi)
			}
		}

		sumRNS[i] = sum
	}

	var smallNorm = true
	for i := 1; i < level+1; i++ {
		smallNorm = smallNorm && (sumRNS[0] == sumRNS[i])
	}

	if !smallNorm {
		var crtReconstruction *big.Int

		sumBigInt := ring.NewUint(0)
		QiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := ringQ.ModulusAtLevel[level]

		for i := 0; i < level+1; i++ {
			QiB.SetUint64(ringQ.Modulus[i])
			crtReconstruction = new(big.Int)
			crtReconstruction.Quo(modulusBigint, QiB)
			tmp.ModInverse(crtReconstruction, QiB)
			tmp.Mod(tmp, QiB)
			crtReconstruction.Mul(crtReconstruction, tmp)

			sumBigInt.Add(sumBigInt, tmp.Mul(ring.NewUint(sumRNS[i]), crtReconstruction))
		}

		sumBigInt.Mod(sumBigInt, modulusBigint)

		logSum = sumBigInt.BitLen()
	} else {
		logSum = bits.Len64(sumRNS[0])
	}

	return
}

func testGenKeyPair(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	// Checks that the secret-key has exactly params.h non-zero coefficients
	t.Run(testString(params, "SK"), func(t *testing.T) {

		skInvNTT := NewSecretKey(params)

		if params.PCount() > 0 {
			params.RingP().InvNTTLvl(sk.Value.P.Level(), sk.Value.P, skInvNTT.Value.P)
			for i := range skInvNTT.Value.P.Coeffs {
				var zeros int
				for j := range skInvNTT.Value.P.Coeffs[i] {
					if skInvNTT.Value.P.Coeffs[i][j] == 0 {
						zeros++
					}
				}
				require.Equal(t, params.ringP.N, zeros+params.h)
			}
		}

		params.RingQ().InvNTTLvl(sk.Value.Q.Level(), sk.Value.Q, skInvNTT.Value.Q)
		for i := range skInvNTT.Value.Q.Coeffs {
			var zeros int
			for j := range skInvNTT.Value.Q.Coeffs[i] {
				if skInvNTT.Value.Q.Coeffs[i][j] == 0 {
					zeros++
				}
			}
			require.Equal(t, params.ringQ.N, zeros+params.h)
		}

	})

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, "PK"), func(t *testing.T) {

		if params.PCount() > 0 {

			params.RingQP().MulCoeffsMontgomeryAndAddLvl(sk.Value.Q.Level(), sk.Value.P.Level(), sk.Value, pk.Value[1], pk.Value[0])
			params.RingQP().InvNTTLvl(sk.Value.Q.Level(), sk.Value.P.Level(), pk.Value[0], pk.Value[0])

			log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].Q.Level(), params.RingQ(), pk.Value[0].Q))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].P.Level(), params.RingP(), pk.Value[0].P))
		} else {
			params.RingQ().MulCoeffsMontgomeryAndAdd(sk.Value.Q, pk.Value[1].Q, pk.Value[0].Q)
			params.RingQ().InvNTT(pk.Value[0].Q, pk.Value[0].Q)

			log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()))
			require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0].Q.Level(), params.RingQ(), pk.Value[0].Q))
		}

	})

}

func testSwitchKeyGen(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	// Checks that switching keys are en encryption under the output key
	// of the RNS decomposition of the input key by
	// 1) Decrypting the RNS decomposed input key
	// 2) Reconstructing the key
	// 3) Checking that the difference with the input key has a small norm
	t.Run(testString(params, "SWKGen/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		ringQ := params.RingQ()
		ringP := params.RingP()
		ringQP := params.RingQP()
		skIn := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		levelQ, levelP := params.QCount()-1, params.PCount()-1

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		swk := NewSwitchingKey(params, levelQ, levelP)
		kgen.(*keyGenerator).genSwitchingKey(skIn.Value.Q, skOut.Value, swk)

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for j := range swk.Value[levelP] {
			ringQP.MulCoeffsMontgomeryAndAddLvl(levelQ, levelP, swk.Value[levelP][j][1], skOut.Value, swk.Value[levelP][j][0])
		}

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [w*P*sOut + e]) = P*sOut + sum(e)
		for j := range swk.Value[levelP] {
			if j > 0 {
				ringQP.AddLvl(levelQ, levelP, swk.Value[levelP][0][0], swk.Value[levelP][j][0], swk.Value[levelP][0][0])
			}
		}

		// sOut * P
		ringQ.MulScalarBigint(skIn.Value.Q, ringP.ModulusAtLevel[levelP], skIn.Value.Q)

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(swk.Value[levelP][0][0].Q, skIn.Value.Q, swk.Value[levelP][0][0].Q)

		// Checks that the error is below the bound
		// Worst error bound is N * floor(6*sigma) * #Keys

		ringQP.InvNTTLvl(levelQ, levelP, swk.Value[levelP][0][0], swk.Value[levelP][0][0])
		ringQP.InvMFormLvl(levelQ, levelP, swk.Value[levelP][0][0], swk.Value[levelP][0][0])

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()*len(swk.Value[levelP])))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(levelQ, ringQ, swk.Value[levelP][0][0].Q))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(levelP, ringP, swk.Value[levelP][0][0].P))

	})
}

func testEncryptor(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	ringQ := params.RingQ()

	t.Run(testString(params, "Encrypt/Pk/MaxLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Pk/MinLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, pk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Sk/MaxLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "Encrypt/Sk/MinLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value.Q, ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, ciphertext.Value[0]))
	})

	t.Run(testString(params, "ShallowCopy/Sk"), func(t *testing.T) {
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.ShallowCopy()
		skEnc1, skEnc2 := enc1.(*skEncryptor), enc2.(*skEncryptor)
		require.True(t, skEnc1.encryptorBase == skEnc2.encryptorBase)
		require.True(t, skEnc1.sk == skEnc2.sk)
		require.False(t, skEnc1.basisextender == skEnc2.basisextender)
		require.False(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.False(t, skEnc1.encryptorSamplers == skEnc2.encryptorSamplers)
	})

	t.Run(testString(params, "ShallowCopy/Pk"), func(t *testing.T) {
		enc1 := NewEncryptor(params, pk)
		enc2 := enc1.ShallowCopy()
		pkEnc1, pkEnc2 := enc1.(*pkEncryptor), enc2.(*pkEncryptor)
		require.True(t, pkEnc1.encryptorBase == pkEnc2.encryptorBase)
		require.True(t, pkEnc1.pk == pkEnc2.pk)
		require.False(t, pkEnc1.basisextender == pkEnc2.basisextender)
		require.False(t, pkEnc1.encryptorBuffers == pkEnc2.encryptorBuffers)
		require.False(t, pkEnc1.encryptorSamplers == pkEnc2.encryptorSamplers)
	})

	sk2 := kgen.GenSecretKey()

	t.Run(testString(params, "WithKey/Sk->Sk"), func(t *testing.T) {
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.WithKey(sk2)
		skEnc1, skEnc2 := enc1.(*skEncryptor), enc2.(*skEncryptor)
		require.True(t, skEnc1.encryptorBase == skEnc2.encryptorBase)
		require.False(t, skEnc1.sk == skEnc2.sk)
		require.False(t, skEnc1.basisextender == skEnc2.basisextender)
		require.False(t, skEnc1.encryptorBuffers == skEnc2.encryptorBuffers)
		require.False(t, skEnc1.encryptorSamplers == skEnc2.encryptorSamplers)
	})

	t.Run(testString(params, "WithKey/Sk->Pk"), func(t *testing.T) {
		enc1 := NewEncryptor(params, sk)
		enc2 := enc1.WithKey(pk)
		skEnc1, pkEnc2 := enc1.(*skEncryptor), enc2.(*pkEncryptor)
		require.True(t, skEnc1.encryptorBase == pkEnc2.encryptorBase)
		require.False(t, skEnc1.basisextender == pkEnc2.basisextender)
		require.False(t, skEnc1.encryptorBuffers == pkEnc2.encryptorBuffers)
		require.False(t, skEnc1.encryptorSamplers == pkEnc2.encryptorSamplers)
	})
}

func testDecryptor(kgen KeyGenerator, t *testing.T) {
	params := kgen.(*keyGenerator).params
	sk := kgen.GenSecretKey()
	ringQ := params.RingQ()
	encryptor := NewEncryptor(params, sk)
	decryptor := NewDecryptor(params, sk)

	t.Run(testString(params, "Decrypt/MaxLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})

	t.Run(testString(params, "Encrypt/MinLevel"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)
		decryptor.Decrypt(ciphertext, plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, plaintext.Value))
	})
}

func testKeySwitcher(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	t.Run(testString(params, "KeySwitch/"), func(t *testing.T) {

		if params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		sk := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		ks := NewKeySwitcher(params)

		ringQ := params.RingQ()
		ringP := params.RingP()

		levelQ := params.MaxLevel()
		alpha := params.PCount()
		levelP := alpha - 1

		QBig := ring.NewUint(1)
		for i := range ringQ.Modulus[:levelQ+1] {
			QBig.Mul(QBig, ring.NewUint(ringQ.Modulus[i]))
		}

		PBig := ring.NewUint(1)
		for i := range ringP.Modulus[:levelP+1] {
			PBig.Mul(PBig, ring.NewUint(ringP.Modulus[i]))
		}

		plaintext := NewPlaintext(params, levelQ)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptor(params, sk)
		ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
		encryptor.Encrypt(plaintext, ciphertext)

		// Tests that a random polynomial decomposed is equal to its
		// reconstruction mod each RNS
		t.Run(testString(params, "DecomposeNTT/"), func(t *testing.T) {

			c2InvNTT := ringQ.NewPolyLvl(ciphertext.Level())
			ringQ.InvNTT(ciphertext.Value[1], c2InvNTT)

			coeffsBigintHaveQ := make([]*big.Int, ringQ.N)
			coeffsBigintHaveP := make([]*big.Int, ringQ.N)
			coeffsBigintRef := make([]*big.Int, ringQ.N)
			coeffsBigintWant := make([]*big.Int, ringQ.N)

			for i := range coeffsBigintRef {
				coeffsBigintHaveQ[i] = new(big.Int)
				coeffsBigintHaveP[i] = new(big.Int)
				coeffsBigintRef[i] = new(big.Int)
				coeffsBigintWant[i] = new(big.Int)
			}

			ringQ.PolyToBigintCenteredLvl(ciphertext.Level(), c2InvNTT, 1, coeffsBigintRef)

			tmpQ := ringQ.NewPolyLvl(ciphertext.Level())
			tmpP := ringP.NewPolyLvl(levelP)

			for i := 0; i < len(ks.BuffDecompQP); i++ {

				ks.DecomposeSingleNTT(levelQ, levelP, alpha, i, ciphertext.Value[1], c2InvNTT, ks.BuffDecompQP[i].Q, ks.BuffDecompQP[i].P)

				// Compute q_alpha_i in bigInt
				qalphai := ring.NewInt(1)

				for j := 0; j < alpha; j++ {
					idx := i*alpha + j
					if idx > levelQ {
						break
					}
					qalphai.Mul(qalphai, ring.NewUint(ringQ.Modulus[idx]))
				}

				ringQ.ReduceLvl(levelQ, ks.BuffDecompQP[i].Q, ks.BuffDecompQP[i].Q)
				ringP.ReduceLvl(levelP, ks.BuffDecompQP[i].P, ks.BuffDecompQP[i].P)

				ringQ.InvNTTLvl(levelQ, ks.BuffDecompQP[i].Q, tmpQ)
				ringP.InvNTTLvl(levelP, ks.BuffDecompQP[i].P, tmpP)

				ringQ.PolyToBigintCenteredLvl(levelQ, tmpQ, 1, coeffsBigintHaveQ)
				ringP.PolyToBigintCenteredLvl(levelP, tmpP, 1, coeffsBigintHaveP)

				// Checks that Reconstruct(NTT(c2 mod Q)) mod q_alpha_i == Reconstruct(NTT(Decomp(c2 mod Q, q_alpha-i) mod QP))
				for i := range coeffsBigintWant[:1] {

					coeffsBigintWant[i].Mod(coeffsBigintRef[i], qalphai)
					coeffsBigintWant[i].Mod(coeffsBigintWant[i], QBig)
					coeffsBigintHaveQ[i].Mod(coeffsBigintHaveQ[i], QBig)
					require.Equal(t, coeffsBigintHaveQ[i].Cmp(coeffsBigintWant[i]), 0)

					coeffsBigintWant[i].Mod(coeffsBigintRef[i], qalphai)
					coeffsBigintWant[i].Mod(coeffsBigintWant[i], PBig)
					coeffsBigintHaveP[i].Mod(coeffsBigintHaveP[i], PBig)
					require.Equal(t, coeffsBigintHaveP[i].Cmp(coeffsBigintWant[i]), 0)

				}
			}
		})

		// Test that Dec(KS(Enc(ct, sk), skOut), skOut) has a small norm
		t.Run(testString(params, "KeySwitch/Standard/"), func(t *testing.T) {
			swk := kgen.GenSwitchingKey(sk, skOut)
			maxLevel := ciphertext.Level()

			for level := maxLevel; level > 0; level-- {
				ks.SwitchKeysInPlace(level, ciphertext.Value[1], swk, ks.BuffQP[1].Q, ks.BuffQP[2].Q)
				ringQ.AddLvl(level, ciphertext.Value[0], ks.BuffQP[1].Q, ks.BuffQP[1].Q)
				ringQ.MulCoeffsMontgomeryAndAddLvl(level, ks.BuffQP[2].Q, skOut.Value.Q, ks.BuffQP[1].Q)
				ringQ.InvNTTLvl(level, ks.BuffQP[1].Q, ks.BuffQP[1].Q)
				require.GreaterOrEqual(t, 11+params.LogN(), log2OfInnerSum(level, ringQ, ks.BuffQP[1].Q),
					fmt.Sprintf("Level: %v LevelSP: %v", level, ks.LevelSP[level]))
			}
		})

	})
}

func testKeySwitchDimension(kgen KeyGenerator, t *testing.T) {

	paramsLargeDim := kgen.(*keyGenerator).params

	t.Run(testString(paramsLargeDim, "KeySwitchDimension/"), func(t *testing.T) {

		if paramsLargeDim.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		paramsSmallDim, _ := NewParametersFromLiteral(ParametersLiteral{
			LogN:     paramsLargeDim.LogN() - 1,
			Q:        paramsLargeDim.Q()[:1],
			P:        paramsLargeDim.P()[:1],
			Sigma:    DefaultSigma,
			RingType: paramsLargeDim.RingType(),
		})

		t.Run("LargeToSmall/", func(t *testing.T) {

			ringQLargeDim := paramsLargeDim.RingQ()
			ringQSmallDim := paramsSmallDim.RingQ()

			kgenLargeDim := NewKeyGenerator(paramsLargeDim)
			skLargeDim := kgenLargeDim.GenSecretKey()
			kgenSmallDim := NewKeyGenerator(paramsSmallDim)
			skSmallDim := kgenSmallDim.GenSecretKey()

			swk := kgenLargeDim.GenSwitchingKey(skLargeDim, skSmallDim)

			plaintext := NewPlaintext(paramsLargeDim, paramsLargeDim.MaxLevel())
			plaintext.Value.IsNTT = true
			encryptor := NewEncryptor(paramsLargeDim, skLargeDim)
			ctLargeDim := NewCiphertextNTT(paramsLargeDim, 1, plaintext.Level())
			encryptor.Encrypt(plaintext, ctLargeDim)

			ks := NewKeySwitcher(paramsLargeDim)
			ks.SwitchKeysInPlace(paramsSmallDim.MaxLevel(), ctLargeDim.Value[1], swk, ks.BuffQP[1].Q, ks.BuffQP[2].Q)
			ringQLargeDim.AddLvl(paramsSmallDim.MaxLevel(), ctLargeDim.Value[0], ks.BuffQP[1].Q, ctLargeDim.Value[0])
			ring.CopyValues(ks.BuffQP[2].Q, ctLargeDim.Value[1])

			//Extracts Coefficients
			ctSmallDim := NewCiphertextNTT(paramsSmallDim, 1, paramsSmallDim.MaxLevel())

			SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQSmallDim, ringQLargeDim, ctSmallDim)

			// Decrypts with smaller dimension key
			ringQSmallDim.MulCoeffsMontgomeryAndAddLvl(ctSmallDim.Level(), ctSmallDim.Value[1], skSmallDim.Value.Q, ctSmallDim.Value[0])
			ringQSmallDim.InvNTTLvl(ctSmallDim.Level(), ctSmallDim.Value[0], ctSmallDim.Value[0])

			require.GreaterOrEqual(t, 10+paramsSmallDim.LogN(), log2OfInnerSum(ctSmallDim.Level(), ringQSmallDim, ctSmallDim.Value[0]))
		})

		t.Run("SmallToLarge/", func(t *testing.T) {

			ringQLargeDim := paramsLargeDim.RingQ()

			kgenLargeDim := NewKeyGenerator(paramsLargeDim)
			skLargeDim := kgenLargeDim.GenSecretKey()
			kgenSmallDim := NewKeyGenerator(paramsSmallDim)
			skSmallDim := kgenSmallDim.GenSecretKey()

			swk := kgenLargeDim.GenSwitchingKey(skSmallDim, skLargeDim)

			plaintext := NewPlaintext(paramsSmallDim, paramsSmallDim.MaxLevel())
			plaintext.Value.IsNTT = true
			encryptor := NewEncryptor(paramsSmallDim, skSmallDim)
			ctSmallDim := NewCiphertextNTT(paramsSmallDim, 1, plaintext.Level())
			encryptor.Encrypt(plaintext, ctSmallDim)

			//Extracts Coefficients
			ctLargeDim := NewCiphertextNTT(paramsLargeDim, 1, plaintext.Level())

			SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, nil, ctLargeDim)

			ks := NewKeySwitcher(paramsLargeDim)
			ks.SwitchKeysInPlace(ctLargeDim.Value[1].Level(), ctLargeDim.Value[1], swk, ks.BuffQP[1].Q, ks.BuffQP[2].Q)
			ringQLargeDim.Add(ctLargeDim.Value[0], ks.BuffQP[1].Q, ctLargeDim.Value[0])
			ring.CopyValues(ks.BuffQP[2].Q, ctLargeDim.Value[1])

			// Decrypts with smaller dimension key
			ringQLargeDim.MulCoeffsMontgomeryAndAddLvl(ctLargeDim.Level(), ctLargeDim.Value[1], skLargeDim.Value.Q, ctLargeDim.Value[0])
			ringQLargeDim.InvNTTLvl(ctLargeDim.Level(), ctLargeDim.Value[0], ctLargeDim.Value[0])

			require.GreaterOrEqual(t, 10+paramsSmallDim.LogN(), log2OfInnerSum(ctLargeDim.Level(), ringQLargeDim, ctLargeDim.Value[0]))
		})
	})
}
