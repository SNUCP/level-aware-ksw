package bench_test

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

func BenchmarkKeySwitch(b *testing.B) {

	level := 27
	pCount := 1

	paramsLitetal := ckks.PN16QP1761
	paramsLitetal.P = paramsLitetal.P[:pCount]

	params, _ := ckks.NewParametersFromLiteral(paramsLitetal)
	prng, _ := utils.NewPRNG()

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKey()
	skOut := kgen.GenSecretKey()
	swk := kgen.GenSwitchingKey(sk, skOut)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{})

	// adjust spindex here
	ksw := evaluator.GetKeySwitcher()
	ksw.SPIndex[level] = 3

	swk = ksw.PreprocessSwitchKey(ksw.SPIndex[level], swk)

	ctIn := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())
	ctOut := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())

	b.Run(fmt.Sprintf("level-%v-%v", level, ksw.LevelPk(level)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.SwitchKeys(ctIn, swk, ctOut)
		}
	})
}
