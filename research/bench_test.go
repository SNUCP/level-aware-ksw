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
		0xfffff4c0001, 0xfffff3a0001, 0xfffff300001, 0xffffee00001,
		0xffffea60001, 0xffffe940001, 0xffffe920001, 0xffffe7c0001,
		0xffffe520001, 0xffffe340001, 0xffffe260001, 0xffffdfc0001,
		0xffffdda0001, 0xffffdd40001, 0xffffdc60001, 0xffffd960001,
		0xffffd800001, 0xffffcc60001, 0xffffc9c0001, 0xffffc780001,
		0xffffc400001, 0xffffc3c0001, 0xffffc2e0001, 0xffffc240001,
		0xffffbc80001, 0xffffbbe0001, 0xffffba00001, 0xffffb7a0001,
		0xffffb5c0001, 0xffffb4c0001, 0xffffb220001, 0xffffb1a0001,
		0xffffade0001, 0xffffada0001, 0xffffaaa0001, 0xffffa8c0001,
		0xffffa600001, 0xffff9c60001, 0xffff9ac0001, 0xffff98a0001,
	},
	P: []uint64{
		// 44 x 4
		0xfffffc60001, 0xfffffac0001, 0xfffff960001, 0xfffff880001,
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
