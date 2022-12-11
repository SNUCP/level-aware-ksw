package bench_test

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var testLevels = []int{30, 28} // L
var pCounts = []int{2, 4, 8}   // k

func BenchmarkKeySwitchPreProcess(b *testing.B) {
	for _, maxLevel := range testLevels {
		for _, pCount := range pCounts {
			// Generate SPIndexes, Choose the fastest one later
			SPIndexes := make([]int, 0)
			for _, v := range []int{1, 2, 4, 8} {
				index := v/pCount - 1
				if (index+1)*pCount != v {
					continue
				}
				SPIndexes = append(SPIndexes, index)
			}

			for _, level := range []int{maxLevel, 24, 16, 8} {
				for _, SPIndex := range SPIndexes {
					testName := fmt.Sprintf("BenchmarkKeySwitchPreProcess-L-%v/l-%v/k-%v/SPIndex-%d", maxLevel, level, pCount, SPIndex)

					paramsLiteral := ckks.PN16QP1761
					paramsLiteral.P = paramsLiteral.P[:pCount]

					params, _ := ckks.NewParametersFromLiteral(paramsLiteral)
					prng, _ := utils.NewKeyedPRNG([]byte{'b', 'y', 't', 'e'})

					kgen := ckks.NewKeyGenerator(params)
					sk := kgen.GenSecretKey()
					skOut := kgen.GenSecretKey()
					swk := kgen.GenSwitchingKey(sk, skOut)
					evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{})

					ksw := evaluator.GetKeySwitcher()
					ksw.SPIndex[level] = SPIndex

					swk = ksw.PreprocessSwitchKey(ksw.SPIndex[level], swk)

					ctIn := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())
					ctOut := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())

					b.Run(testName, func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							evaluator.SwitchKeys(ctIn, swk, ctOut)
						}
					})
				}
			}
		}
	}
}
