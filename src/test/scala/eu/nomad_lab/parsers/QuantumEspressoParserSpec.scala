/*
   Copyright 2016-2017 The NOMAD Developers Group

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */
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
