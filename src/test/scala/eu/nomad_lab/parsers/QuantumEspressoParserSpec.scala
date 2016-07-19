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
    "test1 (pw.x 5.2.1) with json-events" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in", "json-events") must_== ParseResult.ParseSuccess
    }
    "test1 (pw.x 5.2.1) with json" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=dft1.in", "json") must_== ParseResult.ParseSuccess
    }
    "test2 (pw.x 4.0) with json-events" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/scf-ncpp.ref.1193929391", "json-events") must_== ParseResult.ParseSuccess
    }
    "test2 (pw.x 4.0) with json" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/scf-ncpp.ref.1193929391", "json") must_== ParseResult.ParseSuccess
    }
    "test3 (pw.x 5.2.1, binary, USPP) with json-events" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=uspp1.in.1451158822", "json-events") must_== ParseResult.ParseSuccess
    }
    "test3 (pw.x 5.2.1, binary, USPP) with json" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantum-espresso/test/examples/PWSCF/benchmark.out.r11920.inp=uspp1.in.1451158822", "json") must_== ParseResult.ParseSuccess
    }
  }
}
