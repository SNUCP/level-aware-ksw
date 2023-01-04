package bench_test

import (
	"fmt"
	"runtime"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

var qpCount = 40
var pCount = 4

var levels = []int{
	4 - 1, 8 - 1, 12 - 1, 16 - 1,
	20 - 1, 24 - 1, 28 - 1, 32 - 1,
	36 - 1, 37 - 1, 38 - 1, 39 - 1,
}

var paramsLiteral = ckks.ParametersLiteral{
	LogN:     16,
	LogSlots: 15,
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
		0xfffb5e0001, 0xfffb340001, 0xfffb1a0001, 0xfffb0e0001,
	},
	P: []uint64{
		// 44 x 4
		0xffff8a0001, 0xffff820001, 0xffff780001, 0xffff580001,
	},
	DefaultScale: 1 << 44,
	Sigma:        rlwe.DefaultSigma,
	RingType:     ring.Standard,
}

func BenchmarkKeySwitch(b *testing.B) {

	qCount := qpCount - pCount
	paramsLiteral.Q = paramsLiteral.Q[:qCount]
	paramsLiteral.P = paramsLiteral.P[:pCount]
	params, _ := ckks.NewParametersFromLiteral(paramsLiteral)

	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()
	skOut := kgen.GenSecretKey()
	swk := kgen.GenSwitchingKey(sk, skOut)

	encryptor := ckks.NewEncryptor(params, pk)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{})
	ksw := evaluator.GetKeySwitcher()

	for _, level := range levels {
		if level > params.MaxLevel() {
			continue
		}

		testName := fmt.Sprintf("L-%v/l-%v/SP-%v", qCount, level+1, ksw.LevelSP[level]+1)

		ptxt := ckks.NewPlaintext(params, level, params.DefaultScale())
		ctIn := encryptor.EncryptNew(ptxt)
		ctOut := ckks.NewCiphertext(params, 1, level, params.DefaultScale())

		b.Run(testName, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.SwitchKeys(ctIn, swk, ctOut)
			}
		})
	}

	runtime.GC()
}

func BenchmarkCoeffToSlot(b *testing.B) {

	qCount := qpCount - pCount
	paramsLiteral.Q = paramsLiteral.Q[:qCount]
	paramsLiteral.P = paramsLiteral.P[:pCount]
	params, _ := ckks.NewParametersFromLiteral(paramsLiteral)

	CoeffsToSlotsParametersLiteral := advanced.EncodingMatrixLiteral{
		LogN:                params.LogN(),
		LogSlots:            params.LogSlots(),
		Scaling:             1.0 / float64(2*params.Slots()),
		LinearTransformType: advanced.CoeffsToSlots,
		LevelStart:          params.MaxLevel(),
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{

			{params.QiFloat64(params.MaxLevel() - 3)},
			{params.QiFloat64(params.MaxLevel() - 2)},
			{params.QiFloat64(params.MaxLevel() - 1)},
			{params.QiFloat64(params.MaxLevel() - 0)},
		},
	}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, sk)

	// Generates the encoding matrices
	CoeffsToSlotMatrices := advanced.NewHomomorphicEncodingMatrixFromLiteral(CoeffsToSlotsParametersLiteral, encoder)

	// Gets the rotations indexes for CoeffsToSlots
	rotations := CoeffsToSlotsParametersLiteral.Rotations(params.LogN(), params.LogSlots())

	// Generates the rotation keys
	rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)

	// Creates an evaluator with the rotation keys
	eval := advanced.NewEvaluator(params, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

	testName := fmt.Sprintf("BenchmarkCoeffToSlots: qCount:%v pCount:%v", params.QCount(), params.PCount())

	ptxt := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
	ctIn := encryptor.EncryptNew(ptxt)

	b.Run(testName, func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.CoeffsToSlotsNew(ctIn, CoeffsToSlotMatrices)
		}
	})

	runtime.GC()
}

func BenchmarkEvalMod(b *testing.B) {

	qCount := qpCount - pCount
	paramsLiteral.Q = paramsLiteral.Q[:qCount]
	paramsLiteral.P = paramsLiteral.P[:pCount]
	params, _ := ckks.NewParametersFromLiteral(paramsLiteral)

	evm := advanced.EvalModLiteral{
		Q:             0x4000000120001,
		LevelStart:    params.MaxLevel() - 4,
		SineType:      advanced.Cos2,
		MessageRatio:  256.0,
		K:             325,
		SineDeg:       255,
		DoubleAngle:   4,
		ArcSineDeg:    0,
		ScalingFactor: 1 << 44,
	}

	EvalModPoly := advanced.NewEvalModPolyFromLiteral(evm)

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	rlk := kgen.GenRelinearizationKey(sk, 2)
	encryptor := ckks.NewEncryptor(params, sk)
	eval := advanced.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})

	testName := fmt.Sprintf("BenchmarkEvalMod: qCount:%v pCount:%v", params.QCount(), params.PCount())

	ptxt := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
	ctIn := encryptor.EncryptNew(ptxt)

	b.Run(testName, func(b *testing.B) {
		for i := 0; i < b.N; i++ {

			eval.EvalModNew(ctIn, EvalModPoly)
		}
	})

	runtime.GC()
}

func BenchmarkSlotToCoeff(b *testing.B) {

	qCount := qpCount - pCount
	paramsLiteral.Q = paramsLiteral.Q[:qCount]
	paramsLiteral.P = paramsLiteral.P[:pCount]
	params, _ := ckks.NewParametersFromLiteral(paramsLiteral)

	SlotsToCoeffsParametersLiteral := advanced.EncodingMatrixLiteral{
		LogN:                params.LogN(),
		LogSlots:            params.LogSlots(),
		Scaling:             1.0,
		LinearTransformType: advanced.SlotsToCoeffs,
		LevelStart:          params.MaxLevel() - 16,
		BSGSRatio:           2.0,
		BitReversed:         false,
		ScalingFactor: [][]float64{
			{params.QiFloat64(params.MaxLevel() - 19)},
			{params.QiFloat64(params.MaxLevel() - 18)},
			{params.QiFloat64(params.MaxLevel() - 17)},
			{params.QiFloat64(params.MaxLevel() - 16)},
		},
	}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, sk)

	// Generates the encoding matrices
	SlotsToCoeffsMatrix := advanced.NewHomomorphicEncodingMatrixFromLiteral(SlotsToCoeffsParametersLiteral, encoder)

	// Gets the rotations indexes for SlotsToCoeffs
	rotations := SlotsToCoeffsParametersLiteral.Rotations(params.LogN(), params.LogSlots())

	// Generates the rotation keys
	rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)

	// Creates an evaluator with the rotation keys
	eval := advanced.NewEvaluator(params, rlwe.EvaluationKey{Rlk: nil, Rtks: rotKey})

	testName := fmt.Sprintf("BenchmarkSlotToCoeffs: qCount:%v pCount:%v", params.QCount(), params.PCount())

	ptxt := ckks.NewPlaintext(params, params.MaxLevel(), params.DefaultScale())
	ctReal := encryptor.EncryptNew(ptxt)
	ctImag := encryptor.EncryptNew(ptxt)

	b.Run(testName, func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.SlotsToCoeffsNew(ctReal, ctImag, SlotsToCoeffsMatrix)
		}
	})

	runtime.GC()
}
