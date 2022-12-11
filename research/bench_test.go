package bench_test

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var qpCount = 40
var pCounts = []int{1, 2, 4} // k
var levels = []int{31, 23, 15, 7}

func BenchmarkKeySwitchPreProcess(b *testing.B) {
	for _, pCount := range pCounts {
		qCount := qpCount - pCount
		// Generate SPIndexes, Choose the fastest one later
		spIndexes := make([]int, 0)
		for _, v := range []int{1, 2, 4, 8} {
			if v%pCount == 0 {
				spIndexes = append(spIndexes, v/pCount-1)
			}
		}

		paramsLiteral := ckks.PN16QP1761

		paramsLiteral.Q = paramsLiteral.Q[:qCount]
		paramsLiteral.P = paramsLiteral.P[:pCount]

		params, _ := ckks.NewParametersFromLiteral(paramsLiteral)
		prng, _ := utils.NewKeyedPRNG([]byte{'b', 'y', 't', 'e'})

		kgen := ckks.NewKeyGenerator(params)
		sk := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()
		swk := kgen.GenSwitchingKey(sk, skOut)
		evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{})

		ksw := evaluator.GetKeySwitcher()

		for _, level := range levels {
			for _, sp := range spIndexes {
				if level+sp*pCount+1 > qCount {
					continue
				}
				testName := fmt.Sprintf("L-%v/l-%v/k-%v/m-%d", qCount, level+1, pCount, sp+1)

				ksw.SPIndex[level] = sp

				la_swk := ksw.PreprocessSwitchKey(ksw.SPIndex[level], swk)

				ctIn := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())
				ctOut := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())

				b.Run(testName, func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						evaluator.SwitchKeys(ctIn, la_swk, ctOut)
					}
				})
			}
		}
	}
}
