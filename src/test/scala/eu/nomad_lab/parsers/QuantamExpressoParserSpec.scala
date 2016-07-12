package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object QuantamExpressoParserSpec extends Specification {
  "QuantamExpressoParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantam-expresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantam-expresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in", "json") must_== ParseResult.ParseSuccess
    }
  }
}
