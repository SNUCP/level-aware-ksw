package bench_test

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var ciphertextSwitched *ckks.Ciphertext // For preventing optimizations

func BenchmarkKeySwitch(b *testing.B) {
	params, _ := ckks.NewParametersFromLiteral(ckks.PN16QP1761)
	prng, _ := utils.NewPRNG()

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	skOut := kgen.GenSecretKey()
	swk := kgen.GenSwitchingKey(sk, skOut)

	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{})

	for level := params.MaxLevel(); level > 0; level-- {
		b.Run(fmt.Sprintf("level-%v", level), func(b *testing.B) {
			ciphertext := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())
			ciphertextSwitched = ciphertext.CopyNew()
			b.ResetTimer()

			for i := 0; i < b.N; i++ {
				evaluator.SwitchKeys(ciphertext, swk, ciphertextSwitched)
			}
		})
	}
}
