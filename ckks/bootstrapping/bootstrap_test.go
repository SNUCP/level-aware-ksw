package bootstrapping

import (
	"flag"
	"fmt"
	"runtime"
	"sync"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/advanced"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var minPrec float64 = 12.0

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var testBootstrapping = flag.Bool("test-bootstrapping", false, "run the bootstrapping tests (memory intensive)")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

func ParamsToString(params ckks.Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.Beta())
}

func TestBootstrap(t *testing.T) {

	paramSet := defaultParametersLiteral{
		ckks.ParametersLiteral{
			LogN:         16,
			LogSlots:     15,
			DefaultScale: 1 << 40,
			Sigma:        rlwe.DefaultSigma,
			H:            32768,
			Q: []uint64{
				//0x10000000006e0001, // 60 Q0

				0x4000000120001, // 50 Q0

				0x10000140001, // 40
				0xffffe80001,  // 40
				0xffffc40001,  // 40
				0x100003e0001, // 40

				/*
					0xffffb20001,  // 40
					0x10000500001, // 40
					0xffff940001,  // 40
					0xffff8a0001,  // 40

					0xffff820001,  // 40
					0xffff780001,  // 40
					0x10000960001, // 40
					0x10000a40001, // 40
				*/

				// level: 12

				0x7fffe60001, // 39 StC
				0x7fffe40001, // 39 StC
				0x7fffe00001, // 39 StC

				// level: 15

				0x4000000420001,
				0x4000000660001,
				0x40000007e0001,
				0x4000000800001,
				0x3ffffffd20001,
				0x3ffffffb80001,
				0x3fffffed60001,
				0x3fffffec80001,

				/*
					0xfffffffff840001,  // 60 Sine (double angle)
					0x1000000000860001, // 60 Sine (double angle)
					0xfffffffff6a0001,  // 60 Sine
					0x1000000000980001, // 60 Sine
					0xfffffffff5a0001,  // 60 Sine
					0x1000000000b00001, // 60 Sine
					0x1000000000ce0001, // 60 Sine
					0xfffffffff2a0001,  // 60 Sine
				*/

				// level: 23

				0x100000000060001, // 56 CtS
				0xfffffffff00001,  // 56 CtS
				0xffffffffd80001,  // 56 CtS
				0x1000000002a0001, // 56 CtS

				// level: 27
			},
			P: []uint64{
				0x1fffffffffe00001, // Pi 61
				0x1fffffffffc80001, // Pi 61
				//0x1fffffffffb40001, // Pi 61
				//0x1fffffffff500001, // Pi 61

				//0x1fffffffff420001, // Pi 61
				//0x1fffffffff380001, // Pi 61
			},
		},

		Parameters{
			EphemeralSecretWeight: 32,
			SlotsToCoeffsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.SlotsToCoeffs,
				LevelStart:          7,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x7fffe60001},
					{0x7fffe40001},
					{0x7fffe00001},
				},
			},
			EvalModParameters: advanced.EvalModLiteral{
				//Q:             0x10000000006e0001,
				Q:             0x4000000120001, // 50 Q0
				LevelStart:    15,
				SineType:      advanced.Cos1,
				MessageRatio:  256.0,
				K:             16,
				SineDeg:       30,
				DoubleAngle:   3,
				ArcSineDeg:    0,
				ScalingFactor: 1 << 50,
			},
			CoeffsToSlotsParameters: advanced.EncodingMatrixLiteral{
				LinearTransformType: advanced.CoeffsToSlots,
				LevelStart:          19,
				BSGSRatio:           2.0,
				BitReversed:         false,
				ScalingFactor: [][]float64{
					{0x100000000060001},
					{0xfffffffff00001},
					{0xffffffffd80001},
					{0x1000000002a0001},
				},
			},
		},
	}
	ckksParams := paramSet.SchemeParams
	btpParams := paramSet.BootstrappingParams

	EphemeralSecretWeight := btpParams.EphemeralSecretWeight

	// Test original bootstrpping
	ckksParams.H = EphemeralSecretWeight
	btpParams.EphemeralSecretWeight = 0
	params, _ := ckks.NewParametersFromLiteral(ckksParams)
	testbootstrap(params, true, btpParams, t)
	runtime.GC()
}

func testbootstrap(params ckks.Parameters, original bool, btpParams Parameters, t *testing.T) {

	btpType := "Encapsulation/"

	if original {
		btpType = "Original/"
	}

	t.Run(ParamsToString(params, "Bootstrapping/FullCircuit/"+btpType), func(t *testing.T) {

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()

		encoder := ckks.NewEncoder(params)
		encryptor := ckks.NewEncryptor(params, sk)
		decryptor := ckks.NewDecryptor(params, sk)

		evk := GenEvaluationKeys(btpParams, params, sk)

		btp, err := NewBootstrapper(params, btpParams, evk)
		if err != nil {
			panic(err)
		}

		values := make([]complex128, 1<<params.LogSlots())
		for i := range values {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if 1<<params.LogSlots() > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := ckks.NewPlaintext(params, 0, params.DefaultScale())
		encoder.Encode(values, plaintext, params.LogSlots())

		ciphertexts := make([]*ckks.Ciphertext, 2)
		bootstrappers := make([]*Bootstrapper, 2)
		for i := range ciphertexts {
			ciphertexts[i] = encryptor.EncryptNew(plaintext)
			if i == 0 {
				bootstrappers[i] = btp
			} else {
				bootstrappers[i] = bootstrappers[0].ShallowCopy()
			}
		}

		var wg sync.WaitGroup
		wg.Add(2)
		for i := range ciphertexts {
			go func(index int) {
				ciphertexts[index] = bootstrappers[index].Bootstrapp(ciphertexts[index])
				//btp.SetScale(ciphertexts[index], params.Scale())
				wg.Done()
			}(i)
		}
		wg.Wait()

		for i := range ciphertexts {
			verifyTestVectors(params, encoder, decryptor, values, ciphertexts[i], params.LogSlots(), 0, t)
		}
	})
}

func verifyTestVectors(params ckks.Parameters, encoder ckks.Encoder, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, logSlots int, bound float64, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, logSlots, bound)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, precStats.MeanPrecision.Real, minPrec)
	require.GreaterOrEqual(t, precStats.MeanPrecision.Imag, minPrec)
}
