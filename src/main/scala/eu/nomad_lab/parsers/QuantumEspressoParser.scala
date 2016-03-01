package eu.nomad_lab.parsers

import eu.{nomad_lab=>lab}
import eu.nomad_lab.{JsonUtils, DefaultPythonInterpreter}
import org.{json4s => jn}
import scala.collection.breakOut


object QuantumEspressoParser extends SimpleExternalParserGenerator(
  name = "QuantumEspressoParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("QuantumEspressoParser")) ::
      ("parserId" -> jn.JString("QuantumEspressoParser" + lab.QuantumEspressoVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JString(lab.NomadCoreVersionInfo.version)) ::
          (lab.QuantumEspressoVersionInfo.toMap.map{ case (key, value) =>
            (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """^\s*Program (?<programName>\S+)\s+v\.(?<version>\S+)(?:\s+\(svn\s+rev\.\s+(?<revision>\d+)\s*\))?\s+starts[^\n]+
(?:\s*\n?)*This program is part of the open-source Quantum""".r,
  cmd = Seq(DefaultPythonInterpreter.python2Exe(), "${envDir}/parsers/quantum-espresso/parser/parser-quantum-espresso/parser_quantum_espresso.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-quantum-espresso/parser_quantum_espresso.py",
    "parser-quantum-espresso/setup_paths.py",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/quantum_espresso.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-quantum-espresso" -> "parsers/quantum-espresso/parser/parser-quantum-espresso",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
) {

  val fallbackRe = """\s*Program (?<programName>\S+)\s+v\.(?<version>\S+)(?:\s+\(svn\s+rev\.\s+(?<revision>\d+)\s*\))?\s+starts""".r

  override def isMainFile(filePath: String, bytePrefix: Array[Byte], stringPrefix: Option[String]): Option[ParserMatch] = {
    stringPrefix match {
      case Some(str) =>
        mainFileRe.findFirstMatchIn(str) match {
          case Some(m) =>
            val extraInfo: List[(String, jn.JString)] = m.groupNames.map{ (name: String) =>
              name -> jn.JString(m.group(name)) }(breakOut)
            logger.debug(s"$filePath matches parser $name (extraInfo:${JsonUtils.normalizedStr(jn.JObject(extraInfo))}")
            Some(ParserMatch(mainFileMatchPriority, mainFileMatchWeak, jn.JObject(extraInfo)))
          case None =>
            fallbackRe.findFirstMatchIn(str) match {
              case Some(m) =>
                val extraInfo: List[(String, jn.JString)] = m.groupNames.map{ (name: String) =>
              name -> jn.JString(m.group(name)) }(breakOut)
                logger.debug(s"$filePath might match parser $name (extraInfo:${JsonUtils.normalizedStr(jn.JObject(extraInfo))}")
                Some(ParserMatch(mainFileMatchPriority -1, true, jn.JObject(extraInfo)))
              case None =>
                logger.debug(s"$filePath does *NOT* match parser $name")
                None
            }
        }
      case None =>
        logger.warn(s"parser $name asked about $filePath which has no stringPrefix")
        None
    }

  }
}
