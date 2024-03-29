package bench_test

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
	"github.com/tuneinsight/lattigo/v3/utils"
)

var qpCount = 40
var pCount = 1
var levels = []int{
	39 - 1, 38 - 1, 37 - 1, 36 - 1,
	32 - 1, 28 - 1, 24 - 1, 20 - 1,
	16 - 1, 12 - 1, 8 - 1, 4 - 1,
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

	defaultP := make([]uint64, len(paramsLiteral.P))
	copy(defaultP, paramsLiteral.P)

	qCount := qpCount - pCount
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

		if level > params.MaxLevel() {
			continue
		}

		testName := fmt.Sprintf("L-%v/l-%v/SP-%v", qCount, level+1, ksw.LevelSP[level]+1)
		ctIn := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())
		ctOut := ckks.NewCiphertextRandom(prng, params, 1, level, params.DefaultScale())

		b.Run(testName, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				evaluator.SwitchKeys(ctIn, swk, ctOut)
			}
		})
	}
}
