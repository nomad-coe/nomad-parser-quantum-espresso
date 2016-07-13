package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

/**
 * pw.x output Test files:
 * parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in
 * parsers/quantum-espresso/test/examples/PWSCF/scf-ncpp.ref.1193929391
 *
 */

object QuantumEspressoParserSpec extends Specification {
  "QuantumEspressoParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in", "json") must_== ParseResult.ParseSuccess
    }
  }
}
