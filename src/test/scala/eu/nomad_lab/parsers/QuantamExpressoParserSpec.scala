package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object QuantamExpressoParserSpec extends Specification {
  "QuantamExpressoParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantam-expresso/test/examples", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(QuantumEspressoParser, "parsers/quantam-expresso/test/examples", "json") must_== ParseResult.ParseSuccess
    }
  }
}
