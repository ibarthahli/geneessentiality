import org.nspl._
import org.nspl.saddle._
import org.nspl.awtrenderer._
import fileutils._
import stringsplit._
import org.saddle._
import overrepresentation._

object Helpers {

  def readPathogenicVariants(file: String): Seq[String] =
    openSource(file)(_.getLines.map { line =>
      line.split1('\t')(3).splitM('_')(0).split1('.')(0)

    }.toVector)

  def r2(v1: Vec[Double], v2: Vec[Double]) = {
    import org.saddle._
    import org.saddle.linalg._
    val s1 = v1.stdev
    val s2 = v2.stdev

    val cov = v1.demeaned vv v2.demeaned * (1d / (v1.length - 1))
    cov / (s1 * s2)
  }

  def readDiseaseGenes(file: String,
                       symbol2ensg: Map[String, String],
                       hgncTable: Map[String, String]): Set[String] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(2)
      .firstCol("Gene_OMIM_ID")
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("disease genes can't find " + x)
        }

        ensg.getOrElse(x)
      })
      .index
      .toSeq
      .toSet
  }

  def readDiseaseGenesWithCondition(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String]): Series[String, String] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(2)
      .firstCol("OMIM_conditionName")
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("disease genes can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readDiseaseGenesDominant(file: String,
                               symbol2ensg: Map[String, String],
                               hgncTable: Map[String, String]): Set[String] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(2)
      .rfilter(_.get("Inheritance").get == "Autosomal dominant")
      .firstCol("Gene_OMIM_ID")
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("disease genes can't find " + x)
        }

        ensg.getOrElse(x)
      })
      .index
      .toSeq
      .toSet
  }

// def readPathogenicVariants(file:String) : Seq[()]

  def readCCDS(file: String) =
    openSource(file)(_.getLines.map { line =>
      val spl = line.split1('\t')
      // name -> id
      (spl(2), spl(4).split("\\.")(0))
    }.toList)

  def readEnsembleFile(file: String) =
    openSource(file)(
      _.getLines
        .drop(1)
        .map { line =>
          val spl = line.split1('\t')
          (spl(0), spl(1), spl(2))
        }
        .toVector)

  def readHGNCTable(file: String) =
    openSource(file)(_.getLines.flatMap { line =>
      val spl = line.split1('\t')
      val symbollist = {
        val raw = spl(8)
        if (raw.startsWith("\""))
          raw.drop(1).dropRight(1).split("\\|").toList
        else List(raw)
      }

      val prevnamelist = {
        val raw = spl(10)
        if (raw.startsWith("\""))
          raw.drop(1).dropRight(1).split("\\|").toList
        else List(raw)
      }

      val names = spl(1) :: (symbollist ::: prevnamelist)
      val ensg = spl(19)
      names.map(x => x -> ensg)
    }.toMap)

  def readHGNCReverse(file: String) =
    openSource(file)(_.getLines.map { line =>
      val spl = line.split1('\t')
      val symbol = spl(1)
      val ensg = spl(19)
      ensg -> symbol
    }.toMap)

  def readRVIS(file: String,
               symbol2ensg: Map[String, String],
               hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(params = CsvParams(separChar = '\t'))(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(4)
      .firstCol("%RVIS_ExAC_0.05%(AnyPopn)")
      .mapValues(x => 100d - CsvParser.parseDouble(x))
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("rvis can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readRVISUnknown(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String],
      ensgDescription: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    val v1 = CsvParser
      .parse(params = CsvParams(separChar = '\t'))(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(4)
      .firstCol("%RVIS_ExAC_0.05%(AnyPopn)")
      .mapValues(x => 100d - CsvParser.parseDouble(x))

    val v2 = v1.rank() / v1.rank().max.get

    Series(v2, v1.index).filterIx(x => {

      val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

      ensg.isEmpty || ensgDescription
        .get(ensg.get)
        .map(_.isEmpty)
        .getOrElse(true)
    })
  }

  def readMisZ(file: String,
               symbol2ensg: Map[String, String],
               hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("mis_z")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("misz can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readMisZUnknown(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String],
      ensgDescription: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    {
      val v1 = CsvParser
        .parse(CsvFile(file))
        .withColIndex(0)
        .withRowIndex(1)
        .firstCol("mis_z")
        .mapValues(CsvParser.parseDouble)
      val v2 = v1.rank() / v1.rank().max.get

      Series(v2, v1.index).filterIx(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        ensg.isEmpty || ensgDescription
          .get(ensg.get)
          .map(_.isEmpty)
          .getOrElse(true)
      })
    }
  }

  def readSynCount(file: String,
                   symbol2ensg: Map[String, String],
                   hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("n_syn")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("nsyn can't find " + x)
        }

        ensg.getOrElse(x)
      })

  }

  def readDickinsonEG(file: String,
                      symbol2ensg: Map[String, String],
                      hgncTable: Map[String, String]) = {
    import org.saddle.io._
    val s = CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .firstCol("EG")
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("dickinson can't find " + x)
        }

        ensg.getOrElse(x)
      })

    val lethal =
      s.filter(x => x == "Y").index.toSeq.toSet

    val viable = s.filter(_ == "N").index.toSeq.toSet
    (lethal, viable)
  }

  def readPLI(file: String,
              symbol2ensg: Map[String, String],
              hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("pLI")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("pli can't find " + x)
        }

        ensg.getOrElse(x)
      })

  }

  def readPLIUnknown(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String],
      ensgDescription: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    {
      val v1 = CsvParser
        .parse(CsvFile(file))
        .withColIndex(0)
        .withRowIndex(1)
        .firstCol("pLI")
        .mapValues(CsvParser.parseDouble)
      val v2 = v1.rank() / v1.rank().max.get

      Series(v2, v1.index).filterIx(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        ensg.isEmpty || ensgDescription
          .get(ensg.get)
          .map(_.isEmpty)
          .getOrElse(true)
      })
    }
  }

  def readExacNTrunc(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("n_lof")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("pli can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readPRec(file: String,
               symbol2ensg: Map[String, String],
               hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("pRec")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("prec can't find " + x)
        }

        ensg.getOrElse(x)
      })

  }

  def readLofTool(file: String,
                  symbol2ensg: Map[String, String],
                  hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .firstCol("ExACtool percentile")
      .mapValues(x => 1d - CsvParser.parseDouble(x))
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("loftool can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readLofToolUnknown(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String],
      ensgDescription: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    val v1 = CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .firstCol("ExACtool percentile")
      .mapValues(x => 1d - CsvParser.parseDouble(x))

    val v2 = v1.rank() / v1.rank().max.get

    Series(v2, v1.index).filterIx(x => {

      val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

      ensg.isEmpty || ensgDescription
        .get(ensg.get)
        .map(_.isEmpty)
        .getOrElse(true)
    })
  }

  def readPhiPower(file: String,
                   symbol2ensg: Map[String, String],
                   hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(params = CsvParams(separChar = '\t'))(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .firstCol("power")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("phi power can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readPhiUnknown(
      file: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String],
      ensgDescription: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    {
      val v1 = CsvParser
        .parse(params = CsvParams(separChar = '\t'))(CsvFile(file))
        .withColIndex(0)
        .withRowIndex(0)
        .firstCol("phi")
        .mapValues(CsvParser.parseDouble)

      val v2 = v1.rank() / v1.rank().max.get

      Series(v2, v1.index).filterIx(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        ensg.isEmpty || ensgDescription
          .get(ensg.get)
          .map(_.isEmpty)
          .getOrElse(true)
      })
    }
  }

  def readPhi(file: String,
              symbol2ensg: Map[String, String],
              hgncTable: Map[String, String]): Series[String, Double] = {
    import org.saddle.io._
    CsvParser
      .parse(params = CsvParams(separChar = '\t'))(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .firstCol("phi")
      .mapValues(CsvParser.parseDouble)
      .mapIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("phi power can't find " + x)
        }

        ensg.getOrElse(x)
      })
  }

  def readLekMgiEssential(file: String,
                          symbol2ensg: Map[String, String],
                          hgncTable: Map[String, String]) = {
    openSource(file)(_.getLines.toVector.map { x =>
      val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

      if (ensg.isEmpty) {
        println("HART can't find " + x)
      }

      ensg.getOrElse(x)
    })
  }.toSet

  def readHart(file: String,
               symbol2ensg: Map[String, String],
               hgncTable: Map[String, String]) = {
    import org.saddle.io._

    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .col("BF_hct116", "BF_hela", "BF_gbm", "BF_rpe1", "BF_dld1")
      .mapRowIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("HART can't find " + x)
        }

        ensg.getOrElse(x)
      })
      .mapValues(CsvParser.parseDouble)
      .mapColIndex(_ match {
        case "BF_hct116" => "Hart HCT116"
        case "BF_hela" => "Hart HeLa"
        case "BF_gbm" => "Hart GBM"
        case "BF_rpe1" => "Hart RPE1"
        case "BF_dld1" => "Hart DLD1"
      })
  }

  def readHartUnknown(file: String,
                      symbol2ensg: Map[String, String],
                      hgncTable: Map[String, String],
                      ensgDescription: Map[String, String]) = {
    import org.saddle.io._

    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .col("BF_hct116", "BF_hela", "BF_gbm", "BF_rpe1", "BF_dld1")
      .mapValues(CsvParser.parseDouble)
      .mapVec(x => x.rank() / x.rank().max.get)
      .rfilterIx(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        ensg.isEmpty || ensgDescription
          .get(ensg.get)
          .map(_.isEmpty)
          .getOrElse(true)
      })
  }

  def readKeggHierarchy(file: String): Map[String, String] = openSource(file) {
    source =>
      val lines = source.getLines
        .filter(x => x.startsWith("B") || x.startsWith("C"))
        .toVector
      def loop(lines: Vector[String],
               acc: List[(String, Set[String])]): List[(String, Set[String])] =
        if (lines.isEmpty) acc
        else {
          val b = lines.head
          val cs = lines.drop(1).takeWhile(_.startsWith("C"))
          val rest = lines.drop(1 + cs.size)
          loop(
            rest,
            (b.drop(1).trim -> cs
              .map(x =>
                "KEGG_" + x
                  .drop(12)
                  .filterNot(x =>
                    x == '(' || x == ')' || x == '\'' || x == '/' || x == ',')
                  .replaceAllLiterally("-", "_")
                  .replaceAll("\\s+", "_")
                  .replaceAll("_+", "_")
                  .toUpperCase)
              .toSet) :: acc)
        }
      loop(lines, Nil).flatMap(x => x._2.toSeq.map(y => y -> x._1)).toMap
  }

  def readStringTsv(
      f: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String]): Seq[(String, String, Double)] =
    openSource(f)(
      _.getLines
        .drop(1)
        .map { line =>
          val spl = line.split1('\t')
          // database + experimental evidence
          (spl(0), spl(1), math.max(spl(11).toDouble, spl(12).toDouble))
        }
        .toVector).filter(_._3 > 0.7d).map {
      case (s1, s2, v) =>
        (symbol2ensg.get(s1).orElse(hgncTable.get(s1)).getOrElse(s1),
         symbol2ensg.get(s2).orElse(hgncTable.get(s2)).getOrElse(s2),
         v)
    }

  def readWangCS(file: String,
                 symbol2ensg: Map[String, String],
                 hgncTable: Map[String, String]) = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .col("KBM7 CS", "K562 CS", "Jiyoye CS", "Raji CS")
      .mapRowIndex(x => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        if (ensg.isEmpty) {
          println("WANG can't find " + x)
        }

        ensg.getOrElse(x)
      })
      .mapValues(x => -1d * CsvParser.parseDouble(x))
      .mapColIndex(_ match {
        case "KBM7 CS" => "Wang KBM7"
        case "K562 CS" => "Wang K562"
        case "Jiyoye CS" => "Wang Jiyoye"
        case "Raji CS" => "Wang Raji"
      })

  }

  def readWangCSUnknown(file: String,
                        symbol2ensg: Map[String, String],
                        hgncTable: Map[String, String],
                        ensgDescription: Map[String, String]) = {
    import org.saddle.io._
    CsvParser
      .parse(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(0)
      .col("KBM7 CS", "K562 CS", "Jiyoye CS", "Raji CS")
      .mapValues(x => -1d * CsvParser.parseDouble(x))
      .mapVec(x => x.rank() / x.rank().max.get)
      .rfilterIx((x: String) => {

        val ensg = symbol2ensg.get(x).orElse(hgncTable.get(x))

        ensg.isEmpty || ensgDescription
          .get(ensg.get)
          .map(_.isEmpty)
          .getOrElse(true)
      })

  }

  def readBlomenRatio(file: String) = {
    import org.saddle.io._
    CsvParser
      .parse(params = CsvParams(skipLines = 1))(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("ratio")
      .mapValues(x => -1 * (CsvParser.parseDouble(x)))
  }

  def readBlomenRatioUnknown(file: String,
                             ensgDescription: Map[String, String]) = {
    import org.saddle.io._
    CsvParser
      .parse(params = CsvParams(skipLines = 1))(CsvFile(file))
      .withColIndex(0)
      .withRowIndex(1)
      .firstCol("ratio")
      .mapValues(x => -1 * (CsvParser.parseDouble(x)))
      .filterIx(x => ensgDescription.get(x).map(_.isEmpty).getOrElse(true))
  }

  def pearsonsPhi(ct: ContingencyTable2x2): Double = {
    import ct._
    (a11 * a22 - a12 * a21) / scala.math.sqrt(
      rowMargin1 * rowMargin2 * columnMargin2 * columnMargin1)
  }

  def or(ct: ContingencyTable2x2): Double = {
    import ct._
    (a11.toDouble / a12) / (a21.toDouble / a22)
  }

  def readMousePhenotypes(
      impcPhenoGenoFile: String,
      mouse2humanFile: String,
      alleleFile: String,
      symbol2ensg: Map[String, String],
      hgncTable: Map[String, String]
  ) = {
    val mouse2human: Map[String, String] =
      openSource(mouse2humanFile)(_.getLines.map { line =>
        val spl = line.split1('\t')
        spl(5).trim -> spl(0)
      }.toMap)

    val koalleles: Set[String] = openSource(alleleFile)(
      _.getLines
        .filterNot(_.startsWith("#"))
        .map { line =>
          val spl = line.split1('\t')
          spl(0) -> spl(4).split1('|')
        }
        .filter(_._2.contains("Null/knockout"))
        .map(_._1)
        .toSet)

    // marker,allele,zygosity,phenotype
    val impc: Seq[(String, String, String, String)] = {
      import org.saddle.io._
      CsvParser
        .parse(CsvFile(impcPhenoGenoFile))
        .withColIndex(0)
        .col("marker_accession_id",
             "allele_accession_id",
             "zygosity",
             "mp_term_name")
        .toRowSeq
        .map(x => x._2.toVec)
        .map(x => (x.raw(0), x.raw(1), x.raw(2), x.raw(3)))

    }

    impc.map {
      case (marker, allele, zygosity, phenotype) =>
        val humangenesymbol = mouse2human.get(marker)
        val nullallele = koalleles.contains(allele)
        (humangenesymbol, nullallele, zygosity, phenotype)
    }.filter(x => x._2 && x._1.isDefined)
      .map(x => (x._1.get, x._3, x._4))
      .map {
        case (symb, zyg, mp) =>
          val ensg = symbol2ensg.get(symb).orElse(hgncTable.get(symb))

          if (ensg.isEmpty) {
            println("impc can't find " + symb)
          }

          (ensg.getOrElse(symb), zyg, mp)
      }
  }

  def break = throw new RuntimeException

  def readLofs(source: scala.io.Source) =
    source.getLines.dropWhile(_.startsWith("#")).flatMap { line =>
      val columns = line.split1('\t')

      val infoElements = columns(7).split1(';').map { e =>
        val sp = e.split1('=')
        (sp.head, sp.last)
      }

      val pass = columns(6)
      val ref = columns(3)
      val alleles = columns(4).split1(',')

      val allelesAsInVep = alleles.map { alt =>
        if (ref.size == 1 && alt.size == 1) alt
        else if (ref.size == alt.size)
          (ref zip alt)
            .dropWhile(x => x._1 == x._2)
            .reverse
            .dropWhile(x => x._1 == x._2)
            .reverse
            .map(_._2)
            .mkString
        else if (ref.size > alt.size) "-"
        else {
          def loop(r: String, a: String, i: Int): Int =
            if (i >= r.size || i >= a.size) i
            else if (r(i) == a(i)) loop(r, a, i + 1)
            else i

          val i1 = loop(ref, alt, 0)
          val r2 = ref.drop(i1)
          val a2 = alt.drop(i1)
          val i2 = loop(r2.reverse, a2.reverse, 0)
          a2.reverse.drop(i2).reverse
        }

      }

      val csqInfoElem =
        infoElements.find(_._1 == "CSQ").get._2.split1(',').map(_.split1('|'))

      val afElem = infoElements.find(_._1 == "AF").get._2.split1(',')
      val anElem = infoElements.find(_._1 == "AN").get._2

      val nccElem = infoElements.find(_._1 == "NCC").get._2

      val homConsagElem =
        infoElements.find(_._1 == "Hom_CONSANGUINEOUS").get._2.split1(',')

      def extract(s: String) =
        infoElements
          .find(_._1 == s)
          .map(_._2.split1(','))
          .getOrElse(Vector.fill(alleles.size)("0"))

      val hetAfrElem = extract("Het_AFR")
      val hetAmrElem = extract("Het_AMR")
      val hetEasElem = extract("Het_EAS")
      val hetFinElem = extract("Het_FIN")
      val hetNfeElem = extract("Het_NFE")
      val hetOthElem = extract("Het_OTH")
      val hetSasElem = extract("Het_SAS")

      val homAfrElem = extract("Hom_AFR")
      val homAmrElem = extract("Hom_AMR")
      val homEasElem = extract("Hom_EAS")
      val homFinElem = extract("Hom_FIN")
      val homNfeElem = extract("Hom_NFE")
      val homOthElem = extract("Hom_OTH")
      val homSasElem = extract("Hom_SAS")

      if (pass != "PASS") Nil
      else {
        allelesAsInVep.zipWithIndex.flatMap {
          case (alt, altidx) =>
            val csqMatchesAlt = csqInfoElem.filter(_.head == alt)

            if (csqMatchesAlt.isEmpty) {
              println(ref, alleles, allelesAsInVep, csqInfoElem.map(_.head))
            }

            csqMatchesAlt.flatMap { csq =>
              // Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF|context|ancestral

              val ensg = csq(4)
              val lof = csq(56)
              val transcript = csq(6)
              val cpra = columns(0) + "_" + columns(1) + "_" + columns(3) + "_" + alleles(
                  altidx)
              val consequence = csq(1)
              val af = afElem(altidx)
              val homCons = homConsagElem(altidx)
              // val hemiCount =
              //   hemiAfrElem(altidx).toInt +
              //     hemiAmrElem(altidx).toInt +
              //     hemiEasElem(altidx).toInt +
              //     hemiFinElem(altidx).toInt +
              //     hemiNfeElem(altidx).toInt +
              //     hemiOthElem(altidx).toInt +
              //     hemiSasElem(altidx).toInt

              val hetCount =
                hetAfrElem(altidx).toInt +
                  hetAmrElem(altidx).toInt +
                  hetEasElem(altidx).toInt +
                  hetFinElem(altidx).toInt +
                  hetNfeElem(altidx).toInt +
                  hetOthElem(altidx).toInt +
                  hetSasElem(altidx).toInt

              val homCount =
                homAfrElem(altidx).toInt +
                  homAmrElem(altidx).toInt +
                  homEasElem(altidx).toInt +
                  homFinElem(altidx).toInt +
                  homNfeElem(altidx).toInt +
                  homOthElem(altidx).toInt +
                  homSasElem(altidx).toInt

              if (lof == "HC") {
                List(
                  (cpra,
                   ensg,
                   lof,
                   af,
                   homCons,
                   hetCount,
                   homCount,
                   consequence,
                   transcript))
              } else Nil
            }
        }
      }

    }

}
