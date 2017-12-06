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

import eu.{ nomad_lab => lab }
import eu.nomad_lab.{ JsonUtils, DefaultPythonInterpreter }
import org.{ json4s => jn }
import scala.collection.breakOut

object QuantumEspressoParser extends SimpleExternalParserGenerator(
  name = "QuantumEspressoParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("QuantumEspressoParser")) ::
      ("parserId" -> jn.JString("QuantumEspressoParser" + lab.QuantumEspressoVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.QuantumEspressoVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """^\s*Program (?<programName>\S+)\s+v\.(?<version>\S+)(?:\s+\(svn\s+rev\.\s+(?<revision>\d+)\s*\))?\s+starts[^\n]+
(?:\s*\n?)*This program is part of the open-source Quantum""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/quantum-espresso/parser/parser-quantum-espresso/parser_quantum_espresso.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-quantum-espresso/QuantumEspressoCommon.py",
    "parser-quantum-espresso/QuantumEspressoXC.py",
    "parser-quantum-espresso/parser_quantum_espresso.py",
    "parser-quantum-espresso/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/quantum_espresso.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-quantum-espresso" -> "parsers/quantum-espresso/parser/parser-quantum-espresso",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping(),
  metaInfoEnv = Some(lab.meta.KnownMetaInfoEnvs.quantumEspresso)
) {

  val fallbackRe = """\s*Program (?<programName>\S+)\s+v\.(?<version>\S+)(?:\s+\(svn\s+rev\.\s+(?<revision>\d+)\s*\))?\s+starts""".r

  override def isMainFile(filePath: String, bytePrefix: Array[Byte], stringPrefix: Option[String]): Option[ParserMatch] = {
    stringPrefix match {
      case Some(str) =>
        mainFileRe.findFirstMatchIn(str) match {
          case Some(m) =>
            val extraInfo: List[(String, jn.JString)] = m.groupNames.map { (name: String) =>
              name -> jn.JString(m.group(name))
            }(breakOut)
            Some(ParserMatch(mainFileMatchPriority, mainFileMatchWeak, jn.JObject(extraInfo)))
          case None =>
            fallbackRe.findFirstMatchIn(str) match {
              case Some(m) =>
                val extraInfo: List[(String, jn.JString)] = m.groupNames.map { (name: String) =>
                  name -> jn.JString(m.group(name))
                }(breakOut)
                Some(ParserMatch(mainFileMatchPriority - 1, true, jn.JObject(extraInfo)))
              case None =>
                None
            }
        }
      case None =>
        logger.warn(s"parser $name asked about $filePath which has no stringPrefix")
        None
    }

  }
}
