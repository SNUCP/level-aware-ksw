package advanced

import (
	"math"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

func TestHomomorphicMod(t *testing.T) {
	var err error

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping homomorphic mod tests for GOARCH=wasm")
	}

	ParametersLiteral := ckks.ParametersLiteral{
		LogN:         16,
		LogSlots:     15,
		DefaultScale: 1 << 40,
		Sigma:        rlwe.DefaultSigma,
		H:            32768,
		Q: []uint64{
			0x4000000120001, // 50 Q0

			0x10000140001, // 40
			0xffffe80001,  // 40
			0xffffc40001,  // 40

			// SlotToCoeff: depth 4
			0x40001280001,
			0x40001700001,
			0x40001740001,
			0x40001760001,

			// mod eval: depth 12
			0xffffffffffc0001,  // 60 Sine (double angle)
			0x10000000006e0001, // 60 Sine (double angle)
			0xfffffffff840001,  // 60 Sine (double angle)
			0x1000000000860001, // 60 Sine (double angle)
			0xfffffffff6a0001,  // 60 Sine
			0x1000000000980001, // 60 Sine
			0xfffffffff5a0001,  // 60 Sine
			0x1000000000b00001, // 60 Sine
			0x1000000000ce0001, // 60 Sine
			0xfffffffff2a0001,  // 60 Sine
			0xfffffffff240001,  // 60 Sine
			0x1000000000f00001, // 60 Sine

			// CoeffToSlot: depth 4
			0x10000000032a0001,
			0x1000000003360001,
			0x1000000003680001,
			0x1000000003900001,
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
			0x1fffffffffc80001, // Pi 61
			0x1fffffffffb40001, // Pi 61
			0x1fffffffff500001, // Pi 61

			//0x1fffffffff420001, // Pi 61
			//0x1fffffffff380001, // Pi 61
		},
	}
	//testEvalModMarshalling(t)

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	for _, testSet := range []func(params ckks.Parameters, t *testing.T){
		testEvalMod,
	} {
		testSet(params, t)
		runtime.GC()
	}

}

func testEvalModMarshalling(t *testing.T) {
	t.Run("Marshalling", func(t *testing.T) {

		evm := EvalModLiteral{
			Q:             0x80000000080001,
			LevelStart:    12,
			SineType:      Sin,
			MessageRatio:  256.0,
			K:             14,
			SineDeg:       127,
			DoubleAngle:   0,
			ArcSineDeg:    7,
			ScalingFactor: 1 << 55,
		}

		data, err := evm.MarshalBinary()
		assert.Nil(t, err)

		evmNew := new(EvalModLiteral)
		if err := evmNew.UnmarshalBinary(data); err != nil {
			assert.Nil(t, err)
		}
		assert.Equal(t, evm, *evmNew)
	})
}

func testEvalMod(params ckks.Parameters, t *testing.T) {

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	rlk := kgen.GenRelinearizationKey(sk, 2)
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)
	eval := NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})

	/*
		t.Run("SineChebyshevWithArcSine", func(t *testing.T) {

			evm := EvalModLiteral{
				Q:             0x4000000120001, // 50 Q0
				LevelStart:    15,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 60,
			}

			EvalModPoly := NewEvalModPolyFromLiteral(evm)

			values, _, ciphertext := newTestVectorsEvalMod(params, encryptor, encoder, evm, t)

			scale := math.Exp2(math.Round(math.Log2(float64(evm.Q) / evm.MessageRatio)))

			// Scale the message to Delta = Q/MessageRatio
			eval.ScaleUp(ciphertext, math.Round(scale/ciphertext.Scale), ciphertext)

			// Scale the message up to Sine/MessageRatio
			eval.ScaleUp(ciphertext, math.Round((evm.ScalingFactor/evm.MessageRatio)/ciphertext.Scale), ciphertext)

			// Normalization
			eval.MultByConst(ciphertext, 1/(float64(evm.K)*evm.QDiff()), ciphertext)
			eval.Rescale(ciphertext, params.DefaultScale(), ciphertext)

			// EvalMod
			ciphertext = eval.EvalModNew(ciphertext, EvalModPoly)

			// PlaintextCircuit
			//pi2r := 6.283185307179586/complex(math.Exp2(float64(evm.DoubleAngle)), 0)
			for i := range values {
				values[i] -= complex(evm.MessageRatio*evm.QDiff()*math.Round(real(values[i])/(evm.MessageRatio/evm.QDiff())), 0)
				//values[i] = sin2pi2pi(values[i] / complex(evm.MessageRatio*evm.QDiff(), 0)) * complex(evm.MessageRatio*evm.QDiff(), 0) / 6.283185307179586
			}

			verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), 0, t)
		})
	*/

	t.Run("CosOptimizedChebyshevWithArcSine", func(t *testing.T) {

		evm := EvalModLiteral{
			Q:             0x4000000120001,
			LevelStart:    params.MaxLevel() - 4,
			SineType:      Cos2,
			MessageRatio:  256.0,
			K:             325,
			SineDeg:       255,
			DoubleAngle:   4,
			ArcSineDeg:    0,
			ScalingFactor: 1 << 60,
		}

		EvalModPoly := NewEvalModPolyFromLiteral(evm)

		values, _, ciphertext := newTestVectorsEvalMod(params, encryptor, encoder, evm, t)

		scale := math.Exp2(math.Round(math.Log2(float64(evm.Q) / evm.MessageRatio)))

		// Scale the message to Delta = Q/MessageRatio
		eval.ScaleUp(ciphertext, math.Round(scale/ciphertext.Scale), ciphertext)

		// Scale the message up to Sine/MessageRatio
		eval.ScaleUp(ciphertext, math.Round((evm.ScalingFactor/evm.MessageRatio)/ciphertext.Scale), ciphertext)

		// Normalization
		eval.MultByConst(ciphertext, 1/(float64(evm.K)*evm.QDiff()), ciphertext)
		eval.Rescale(ciphertext, params.DefaultScale(), ciphertext)

		// EvalMod
		ciphertext = eval.EvalModNew(ciphertext, EvalModPoly)

		// PlaintextCircuit
		//pi2r := 6.283185307179586/complex(math.Exp2(float64(evm.DoubleAngle)), 0)
		for i := range values {
			values[i] -= complex(evm.MessageRatio*evm.QDiff()*math.Round(real(values[i])/(evm.MessageRatio/evm.QDiff())), 0)
		}

		verifyTestVectors(params, encoder, decryptor, values, ciphertext, params.LogSlots(), 0, t)
	})
}

func newTestVectorsEvalMod(params ckks.Parameters, encryptor ckks.Encryptor, encoder ckks.Encoder, evm EvalModLiteral, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext) {

	logSlots := params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	K := float64(evm.K - 1)
	Q := float64(evm.Q) / math.Exp2(math.Round(math.Log2(float64(evm.Q)))) * evm.MessageRatio

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(math.Round(utils.RandFloat64(-K, K))*Q+utils.RandFloat64(-1, 1), 0)
	}

	values[0] = complex(K*Q+0.5, 0)

	plaintext = ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())

	encoder.Encode(values, plaintext, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}
