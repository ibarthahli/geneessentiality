import org.nspl._
import org.nspl.saddle._
import org.nspl.awtrenderer._
import fileutils._
import stringsplit._
import org.saddle._
import org.saddle.linalg._

import OverRep._

object Runner {
  def run = {
    import Helpers._

    // val phiFile = "exac_phi.txt"
    val originalPhiFile = "newestimation.absolute.posteriors.txt"
    val loftoolFile = "loftool_fadista.csv.txt"
    val exacFile = "nature19057-SI Table 13.csv.txt"
    val rvisFile = "RVIS_Unpublished_ExAC_May2015.txt"
    val hgncFile = "non_alt_loci_set.txt"
    val ensembleFile = "ensg_enst_name_symbol.txt"
    val blomenFile = "aac7557_SM_Table_S2.csv"
    val wangFile = "aac7041_SM_Table_S3.csv"
    val hartFile = "mmc3.csv"
    val keggHierarchyFile = "br08901.keg"
    val dickinsonFile = "nature19356-s7.csv"
    val diseaseGenesFile = "WGS_gene_condition_v2.0.csv"
    val lekMgiEssentialFile = "mgi_essential.tsv"
    val geneDescriptions = "ensg_description.txt"
    val stringCoreTsv = "core_string_interaction.tsv"
    val exacVCF = "ExAC.r0.3.1.sites.vep.vcf.gz"
    val pathogenicVariantsFile =
      "pathogenic_2016_07_R1_SNVs_noBenign_withInDel_uniqPos_insideENSTmatchedFromHg19.bed"

    // openSource(exacVCF) { s =>
    //   openFileWriter(new java.io.File("lofs.exac.tsv")) { writer =>
    //     readLofs(s).foreach { c =>
    //       writer.write(c.productIterator.mkString("\t"))
    //       writer.write("\n")
    //     }
    //   }
    // }

    val keggH = readKeggHierarchy(keggHierarchyFile)

    val hgncTable = readHGNCTable(hgncFile)
    val ensg2symbol =
      readEnsembleFile(ensembleFile).map(x => x._1 -> x._3).toMap

    val unknowngenes =
      readEnsembleFile(geneDescriptions).map(x => x._1 -> x._3).toMap

    val symbol2ensg = ensg2symbol.map(_.swap)

    val symbol2ensgBoth = symbol2ensg ++ Map(
        "CCL4L2" -> "ENSG00000276125",
        "ABP1" -> "ENSG00000002726",
        "GSTT1" -> "ENSG00000277656",
        "LILRA3" -> "ENSG00000273884",
        "SF3B14" -> "ENSG00000115128",
        "NSG1" -> "ENSG00000168824",
        "HMP19" -> "ENSG00000170091",
        "SGK223" -> "ENSG00000275342",
        "SGK494" -> "ENSG00000167524"
      )
    val ensg = readEnsembleFile(ensembleFile).map(_._1).toSet

    val kegg =
      openSource("c2.cp.kegg.v5.1.symbols.gmt")(overrepresentation.readGMT).map {
        case (name, set) =>
          (name, set.flatMap(x => symbol2ensgBoth.get(x)))
      }
    val keggM = kegg.toMap

    val goBP =
      openSource("c5.bp.v5.1.symbols.gmt")(overrepresentation.readGMT).map {
        case (name, set) =>
          (name, set.flatMap(x => symbol2ensgBoth.get(x)))
      }
    val goBPM = kegg.toMap
    val goBPShort = goBP.filter(x => x._2.size < 100 && x._2.size > 10)
    println("GOBPSHORT SIZE: " + goBPShort.size)

    def printOverrep(set: Set[String],
                     pathways: Seq[(String, Set[String])] = kegg,
                     superPathways: Map[String, String] = keggH,
                     pThreshold: Double = 1E-4,
                     bg: Option[Set[String]] = None) = {
      val bg1 = bg.getOrElse(pathways.map(_._2).reduce(_ ++ _))
      val (enr, depl) = overrep(set, pathways, bg1)
      println(
        "Enrichment:\n" +
          enr
            .filter(_._2 < pThreshold)
            .map(x => superPathways.get(x._1).getOrElse("Unknown"))
            .groupBy(x => x)
            .toSeq
            .map(x => x._1 + ":" + x._2.size)
            .mkString("\n") + "\n" +
          enr
            .filter(_._2 < pThreshold)
            .map(y => y._1 + " " + y._2 + " " + y._3.map(ensg2symbol))
            .mkString("\n") + "\nDepletion:\n" +
          depl
            .filter(_._2 < pThreshold)
            .map(x => superPathways.get(x._1).getOrElse("Unknown"))
            .groupBy(x => x)
            .toSeq
            .map(x => x._1 + ":" + x._2.size)
            .mkString("\n") + "\n" +
          depl
            .filter(_._2 < pThreshold)
            .map(y => y._1 + " " + y._2 + " " + y._3.map(ensg2symbol))
            .mkString("\n")
      )
    }

    def dropDup(series: Series[String, Double], name: String) = {

      val dups =
        series.toSeq.groupBy(_._1).filter(_._2.size > 1).toSeq.map(_._1).toSet

      series.toSeq
        .groupBy(_._1)
        .filter(_._2.size == 1)
        .toSeq
        .filter(_._2.size == 1)
        .map(x => x._2.head._1 -> x._2.head._2)
        .toSeries -> dups
    }
    def dropDupF(f: Frame[String, String, Double], name: String) = {
      val dups =
        f.rowIx.toSeq
          .groupBy(x => x)
          .filter(_._2.size > 1)
          .toSeq
          .map(_._1)
          .toSet

      f.rfilterIx(rx => !dups.contains(rx)) -> dups

    }

    println(
      "Number of >0.9 pLI without any id mapping: " + (readPLI(
        exacFile,
        symbol2ensg,
        hgncTable).filter(_ > 0.9).length))

    val (rvis, rvisdup) =
      dropDup(readRVIS(rvisFile, symbol2ensg, hgncTable), "rvis")
    val (misz, miszdup) =
      dropDup(readMisZ(exacFile, symbol2ensg, hgncTable), "misz")
    val (pli, plidup) =
      dropDup(readPLI(exacFile, symbol2ensg, hgncTable), "pli")
    val (phi, phidup) =
      dropDup(readPhi(originalPhiFile, symbol2ensg, hgncTable), "phi")
    val (loftool, loftooldup) =
      dropDup(readLofTool(loftoolFile, symbol2ensg, hgncTable), "loftool")
    val (blomen, blomendup) = dropDup(readBlomenRatio(blomenFile), "blomen")
    val (wang, wangdup) =
      dropDupF(readWangCS(wangFile, symbol2ensg, hgncTable), "wang")
    val (hart, hartdup) =
      dropDupF(readHart(hartFile, symbol2ensg, hgncTable), "hart")
    val (lofcount, lofdup) =
      dropDup(readExacNTrunc(exacFile, symbol2ensg, hgncTable), "lof")

    val lekMgiEssential =
      readLekMgiEssential(lekMgiEssentialFile, symbol2ensg, hgncTable)

    val (dickinsonEssential1, dickinsonNonEssential1) =
      readDickinsonEG(dickinsonFile, symbol2ensg, hgncTable)

    val diseaseGenes =
      readDiseaseGenes(diseaseGenesFile, symbol2ensg, hgncTable)

    val dominantDiseaseGenes =
      readDiseaseGenesDominant(diseaseGenesFile, symbol2ensg, hgncTable)

    val recessiveDiseaseGenes = diseaseGenes &~ dominantDiseaseGenes

    val diseaseGenesWithCondition =
      readDiseaseGenesWithCondition(diseaseGenesFile, symbol2ensg, hgncTable)

    val phiPower = readPhiPower(originalPhiFile, symbol2ensg, hgncTable)

    val prec = readPRec(exacFile, symbol2ensg, hgncTable)

    println(
      "Total drop b/c dup: " + (rvisdup ++ miszdup ++ plidup ++ phidup ++ loftooldup ++ blomendup ++ wangdup ++ hartdup).size)

    val table = (Frame(
        "RVIS" -> rvis,
        "missense-Z" -> misz,
        "pLi" -> pli,
        "Phi" -> phi,
        "LofTool" -> loftool,
        "Blomen KBM7" -> blomen
      ) rconcat wang rconcat hart) //.row(phiPower.filter(_ > 0.2).index.toSeq: _*)

    val dickinsonEssential = dickinsonEssential1
    val dickinsonNonEssential = dickinsonNonEssential1

    {
      import org.saddle.io.CsvImplicits._
      val tableToPrint = table.rconcat(
        Frame("Dickinson_EG" -> Series(
          dickinsonEssential.map(x => x -> 1.0).toSeq ++ dickinsonNonEssential
            .map(x => x -> 0d): _*)))
      tableToPrint.writeCsvFile("table.csv")
    }

    val invivoScoreNames = Set(
      "RVIS",
      "missense-Z",
      "LofTool",
      "pLi",
      "Phi"
    )

    println(
      "DUPS: " + table.rowIx.toSeq
        .groupBy(x => x)
        .toSeq
        .map(x => x._1 -> x._2.size)
        .filter(_._2 > 1))

    println(
      "MISSING: " + table.rowIx.toSeq
        .filterNot(g => ensg.contains(g))
        .distinct
        .sorted
        .size)

    val percentiles = table.mapVec(x => x.rank() / x.rank().max.get)

    val percentileCutoff = 0.85

    val invitroScoreNames = table.colIx.toSeq.toSet &~ invivoScoreNames

    val medianPercentiles = {
      val invivo = Series(
        percentiles.col(invivoScoreNames.toSeq: _*).toRowSeq.map {
          case (rs, series) =>
            rs -> series.toVec.dropNA.median
        }: _*)

      val invitro = Series(
        percentiles.col(invitroScoreNames.toSeq: _*).toRowSeq.map {
          case (rs, series) =>
            rs -> series.toVec.dropNA.median
        }: _*)

      Frame("invivo" -> invivo, "invitro" -> invitro).rdropNA
    }

    {
      val pubmedcounts = openSource("pubmed_gene_count.csv")(
        _.getLines
          .map(_.split(","))
          .map(x => x(0) -> x(1).toInt)
          .toList
          .toSeries)
      val median = pubmedcounts.toVec.median
      val poorlyStudied = pubmedcounts.filter(_ < 1)

      val essentialPoorlyStudied = medianPercentiles
        .row(poorlyStudied.mapIndex(symbol2ensg).index.toSeq: _*)
        .rfilter(_.toVec.toSeq.exists(_ > 0.95))
        .mapRowIndex(ensg2symbol)

      println(
        "Essential and poorly studied genes: \n" +
          essentialPoorlyStudied.rowIx.toSeq
            .map(x => x -> unknowngenes.get(symbol2ensg(x)))
            .mkString("\n"))
      essentialPoorlyStudied.print(nrows = 31, ncols = 2)
    }

    {
      import org.saddle.io.CsvImplicits._
      medianPercentiles.writeCsvFile("medianpercentiles.csv")
    }

    val invivoScores =
      table.colIx.toSeq.filter(x => invivoScoreNames.contains(x))
    val invitroScores =
      table.colIx.toSeq.filterNot(x => invivoScoreNames.contains(x))

    println("Joined table: ")
    println(table)

    println("Fully joined genes: " + table.rfilter(!_.hasNA).numRows)

    println("Genes with median percentiles: " + medianPercentiles.numRows)

    println(medianPercentiles)

    val r2med = r2(medianPercentiles.firstCol("invivo").toVec,
                   medianPercentiles.firstCol("invitro").toVec)
    println("R2 of median percentiles: " + r2med)

    val invivoEssentialGenes = medianPercentiles
      .firstCol("invivo")
      .filter(_ > percentileCutoff)
      .index
      .toSeq
      .toSet
    val invivoNotEssentialGenes = medianPercentiles
      .firstCol("invivo")
      .filter(_ <= percentileCutoff)
      .index
      .toSeq
      .toSet
    val invitroEssentialGenes = medianPercentiles
      .firstCol("invitro")
      .filter(_ > percentileCutoff)
      .index
      .toSeq
      .toSet
    val invitroNotEssentialGenes = medianPercentiles
      .firstCol("invitro")
      .filter(_ <= percentileCutoff)
      .index
      .toSeq
      .toSet

    println("invivo: " + invivoEssentialGenes.size)

    println("not invivo: " + invivoNotEssentialGenes.size)

    println("invitro: " + invitroEssentialGenes.size)

    println("not invitro: " + invitroNotEssentialGenes.size)

    {
      println("PATHOGENIC VARIANTS")
      val genesWithPathogenicVariants: Seq[(String, Int)] =
        readPathogenicVariants(pathogenicVariantsFile)
          .groupBy(x => x)
          .toSeq
          .map(x => x._1 -> x._2.size)

      val essentialGenesAndDiseases = genesWithPathogenicVariants
        .filter(
          x =>
            invivoEssentialGenes.contains(x._1) || invitroEssentialGenes
              .contains(x._1) || dickinsonEssential.contains(x._1))
        .map {
          case (gene, numVar) =>
            (ensg2symbol(gene),
             diseaseGenesWithCondition.get(gene): Option[String],
             numVar,
             invivoEssentialGenes.contains(gene),
             invitroEssentialGenes.contains(gene),
             dickinsonEssential.contains(gene))
        }
        .filter(_._2.isDefined)

      println(
        "Diseases in in vitro and mouse essential genes (not human): " + essentialGenesAndDiseases
          .filter(x => !x._4 && (x._5 || x._6))
          .map(_._2)
          .toSet
          .size)

      writeToFile(
        "essentialgenes_diseases.tsv",
        "symbol\tcondition\tnumber_of_pathogenic_variants\tinvivo_essential\tinvitro_essential\tmice_essential\n" +
          essentialGenesAndDiseases
            .map(x =>
              (List(x._1, x._2.get, x._3, x._4, x._5, x._6).mkString("\t")))
            .mkString("\n"))

    }

    {
      println("COMMON LOF")
      case class LofRecord(cpra: String,
                           gene: String,
                           af: Double,
                           het: Int,
                           hom: Int,
                           effect: String,
                           transcript: String)

      val lofsWithTranscripts: Seq[LofRecord] =
        openSource("lofs.exac.tsv")(_.getLines.map { line =>
          val spl = line.split1('\t')
          LofRecord(spl(0),
                    spl(1),
                    spl(3).toDouble,
                    spl(5).toInt,
                    spl(6).toInt,
                    spl(7),
                    spl(8))
        }.toList)

      val lofsWithGene = lofsWithTranscripts
        .groupBy(x => (x.cpra, x.gene))
        .toSeq
        .map(x => x._2.head)

      val lofs = lofsWithGene.groupBy(_.cpra).toSeq.map(_._2.head)

      val homozygousLofs = lofs.filter(_.hom > 0)

      println("lof variants af<0.001: " + lofs.count(_.af < 0.001))

      println("lof median af all: " + lofs.map(_.af).toVec.median)

      println(
        "hom lof median af all: " + homozygousLofs.map(_.af).toVec.median)

      println("hom lof variants af < 0.001: " + homozygousLofs.count(x =>
        x.af < 0.001))

      println("hom lof variants af > 0.05: " + homozygousLofs.count(x =>
        x.af > 0.05))

      println("hom lof variants af>0.01 < 0.05: " + homozygousLofs.count(x =>
        x.af <= 0.05 && x.af > 0.01))

      println("hom lof variants af>0.001 < 0.01: " + homozygousLofs.count(x =>
        x.af <= 0.01 && x.af > 0.001))

      println("all hom lof variants: " + homozygousLofs.size)

      val genesWithHomozygousLof =
        lofsWithGene.filter(_.hom > 0).map(_.gene).toSet

      println(
        "genes with homozygous lofs ('loffy'): " + genesWithHomozygousLof.size)

      val genesWithoutHomozygousLofs = lofsWithGene
          .map(_.gene)
          .toSet &~ genesWithHomozygousLof

      println(
        "in vivo essential & loffy: " + (genesWithHomozygousLof & invivoEssentialGenes).size)
      println(
        "in vitro essential & loffy: " + (genesWithHomozygousLof & invitroEssentialGenes).size)
      println(
        "mice essential & loffy: " + (genesWithHomozygousLof & dickinsonEssential).size)

      println(
        "mice essential & in vivo essential & loffy: " + (genesWithHomozygousLof & dickinsonEssential & invivoEssentialGenes).size)

      println(
        "mice essential & in vitro essential & loffy: " + (genesWithHomozygousLof & dickinsonEssential & invitroEssentialGenes).size)

      println(
        "in vivo & in vitro essential & loffy: " + (genesWithHomozygousLof & invivoEssentialGenes & invitroEssentialGenes).size)

      println(
        "in vivo & in vitro essential & mice & loffy: " + (genesWithHomozygousLof & invivoEssentialGenes & dickinsonEssential & invitroEssentialGenes).size)

      println(
        "in vivo & in vitro essential & mice & loffy: " + (genesWithHomozygousLof & invivoEssentialGenes & dickinsonEssential & invitroEssentialGenes))

      val coreEssentialAndLoffy = (genesWithHomozygousLof & invivoEssentialGenes & dickinsonEssential & invitroEssentialGenes)

      println(
        "Core essential genes with homozygous lofs: \n " + coreEssentialAndLoffy
          .map(x => ensg2symbol.get(x).getOrElse(x)))

      println(
        "core essential and loffy variants: \n" + lofsWithGene
          .filter(x => coreEssentialAndLoffy.contains(x.gene) && (x.hom > 0))
          .sortBy(_.gene)
          .mkString("\n"))

      println("genes with homozygous lofs")
      println(
        "median pli: " + genesWithHomozygousLof.toSeq
          .flatMap(x => pli.get(x): Option[Double])
          .toVec
          .median + ", n=" + genesWithHomozygousLof.toSeq
          .flatMap(x => pli.get(x): Option[Double])
          .toVec
          .length)

      println(
        "median prec: " + genesWithHomozygousLof.toSeq
          .flatMap(x => prec.get(x): Option[Double])
          .toVec
          .median)

      println(
        "median median in vivo: " + genesWithHomozygousLof.toSeq
          .flatMap(x =>
            medianPercentiles.firstCol("invivo").get(x): Option[Double])
          .toVec
          .median + ", n=" + genesWithHomozygousLof.toSeq
          .flatMap(x =>
            medianPercentiles.firstCol("invivo").get(x): Option[Double])
          .toVec
          .length)

      println(
        "median median in vitro: " + genesWithHomozygousLof.toSeq
          .flatMap(x =>
            medianPercentiles.firstCol("invitro").get(x): Option[Double])
          .toVec
          .median + ", n=" + genesWithHomozygousLof.toSeq
          .flatMap(x =>
            medianPercentiles.firstCol("invitro").get(x): Option[Double])
          .toVec
          .length)

      println("other genes")
      println(
        "median pli: " + genesWithoutHomozygousLofs.toSeq
          .flatMap(x => pli.get(x): Option[Double])
          .toVec
          .median + ", n=" + genesWithoutHomozygousLofs.toSeq
          .flatMap(x => pli.get(x): Option[Double])
          .toVec
          .length)

      println(
        "median prec: " + genesWithoutHomozygousLofs.toSeq
          .flatMap(x => prec.get(x): Option[Double])
          .toVec
          .median)

      println(
        "median median in vivo: " + genesWithoutHomozygousLofs.toSeq
          .flatMap(x =>
            medianPercentiles.firstCol("invivo").get(x): Option[Double])
          .toVec
          .median)

      println(
        "median median in vitro: " + genesWithoutHomozygousLofs.toSeq
          .flatMap(x =>
            medianPercentiles.firstCol("invitro").get(x): Option[Double])
          .toVec
          .median)

      printOverrep(genesWithHomozygousLof)

    }

    {
      println("DICKINSON ESSENTIAL LIST ")
      println("dickinsonEssential: " + dickinsonEssential.size)
      println("invivoEssentialGenes: " + invivoEssentialGenes.size)
      println(
        "dickinsonNonEssential & invivoNonEssentialGenes: " + (dickinsonNonEssential & invivoNotEssentialGenes).size)
      println(
        "dickinsonEssential & invivoEssentialGenes: " + (dickinsonEssential & invivoEssentialGenes).size)
      println(
        "dickinsonEssential &~ invivoEssentialGenes: " + (dickinsonEssential &~ invivoEssentialGenes).size)
      println(
        "invivoEssentialGenes &~ dickinsonEssential: " + (invivoEssentialGenes &~ dickinsonEssential).size)

      println(
        "(dickinsonEssential &~ in vivo essential) & high power: " + ((dickinsonEssential &~ invivoEssentialGenes) & phiPower
          .filter(_ > 0.2)
          .index
          .toSeq
          .toSet).size)

      println("GENE SETS")
      println("dickinsonEssential & invivoEssentialGenes GO ")
      printOverrep(dickinsonEssential & invivoEssentialGenes, goBP)
      println("dickinsonNonEssential & invivoEssentialGenes")
      printOverrep(dickinsonNonEssential & invivoEssentialGenes, kegg)
      println("dickinsonEssential & invivoEssentialGenes")
      printOverrep(dickinsonEssential & invivoEssentialGenes, kegg)
      println("dickinsonEssential &~ invivoEssentialGenes")
      printOverrep(dickinsonEssential &~ invivoEssentialGenes, kegg)
      println("dickinsonEssential &~ invivoEssentialGenes & high power")
      printOverrep((dickinsonEssential &~ invivoEssentialGenes) & phiPower
                     .filter(_ > 0.2)
                     .index
                     .toSeq
                     .toSet,
                   kegg)

      println("dominant disease genes: " + dominantDiseaseGenes.size)
      println("recessive disease genes: " + recessiveDiseaseGenes.size)
      println(
        " genes dickinsonEssential & invivoEssentialGenes & disease: " + (dickinsonEssential & invivoEssentialGenes & diseaseGenes).size)
      println(
        " genes dickinsonEssential & invivoEssentialGenes & disease dominant: " + (dickinsonEssential & invivoEssentialGenes & dominantDiseaseGenes).size)
      println(
        " genes dickinsonEssential & invivoEssentialGenes & disease recessive: " + (dickinsonEssential & invivoEssentialGenes & recessiveDiseaseGenes).size)

      println(" dickinsonEssential & invivoEssentialGenes enrichment in dominant diseases: " + (overrepresentation
        .enrichmentTest(
          total = table.rowIx.toSeq.toSet.size,
          marked = (dominantDiseaseGenes & (table.rowIx.toSeq.toSet)).size,
          draws = (dickinsonEssential & invivoEssentialGenes).size,
          markedDraws =
            (dickinsonEssential & invivoEssentialGenes & dominantDiseaseGenes).size)))

      println(" dickinsonEssential & invivoEssentialGenes depletion in dominant diseases: " + (overrepresentation
        .depletionTest(
          total = table.rowIx.toSeq.toSet.size,
          marked = (dominantDiseaseGenes & (table.rowIx.toSeq.toSet)).size,
          draws = (dickinsonEssential & invivoEssentialGenes).size,
          markedDraws =
            (dickinsonEssential & invivoEssentialGenes & dominantDiseaseGenes).size)))

      println(" dickinsonEssential & invivoEssentialGenes enrichment in recessive diseases: " + (overrepresentation
        .enrichmentTest(
          total = table.rowIx.toSeq.toSet.size,
          marked = (recessiveDiseaseGenes & (table.rowIx.toSeq.toSet)).size,
          draws = (dickinsonEssential & invivoEssentialGenes).size,
          markedDraws =
            (dickinsonEssential & invivoEssentialGenes & recessiveDiseaseGenes).size)))

      val dickinsonOnlyHighPower = (dickinsonEssential &~ invivoEssentialGenes) & phiPower
          .filter(_ > 0.2)
          .index
          .toSeq
          .toSet

      println(" dickinsonEssential &~ invivoEssentialGenes high power enrichment in recessive diseases: " + (overrepresentation
        .enrichmentTest(
          total = table.rowIx.toSeq.toSet.size,
          marked = (recessiveDiseaseGenes & (table.rowIx.toSeq.toSet)).size,
          draws = (dickinsonOnlyHighPower).size,
          markedDraws = (dickinsonOnlyHighPower & recessiveDiseaseGenes).size)))

      println("DICKINSON ESSENTIAL LIST END \n\n")

    }

    println(
      " invivoEssentialGenes enrichment in dominant diseases: " + (overrepresentation
        .enrichmentTest(
          total = table.rowIx.toSeq.toSet.size,
          marked = (dominantDiseaseGenes & (table.rowIx.toSeq.toSet)).size,
          draws = (invivoEssentialGenes).size,
          markedDraws = (invivoEssentialGenes & dominantDiseaseGenes).size)))

    def pathwayplot1(horizontal: Set[String],
                     vertical: Set[String],
                     pw: Seq[(String, Set[String])] = kegg,
                     background: Set[String]) = {
      val invivoOverRep = overrep(horizontal, pw, background)
      val invitroOverRep = overrep(vertical, pw, background)
      val x = (invivoOverRep._1.map(x => x._1 -> (-1 * math.log10(x._2))) ++
        invivoOverRep._2.map(x => x._1 -> (math.log10(x._2))))
        .groupBy(_._1)
        .toSeq
        .map(x => x._1 -> x._2.map(_._2).maxBy(math.abs))
        .toSeries
      val y = (invitroOverRep._1.map(x => x._1 -> (-1 * math.log10(x._2))) ++
        invitroOverRep._2.map(x => x._1 -> (math.log10(x._2))))
        .groupBy(_._1)
        .toSeq
        .map(x => x._1 -> x._2.map(_._2).maxBy(math.abs))
        .toSeries

      val numTest = math.log10(0.05 / (x.length * 4)) * -1
      val frame =
        Frame("x" -> x, "y" -> y)

      val colors = frame
        .rfilter(
          series =>
            math.abs(series.toVec.raw(0)) > numTest || math.abs(
              series.toVec.raw(1)) > numTest)
        .rowIx
        .toSeq

      (frame, colors)
    }

    def pathwayplot2(frame: Frame[String, String, Double],
                     colored: Seq[String],
                     legend: Seq[String],
                     labelsToDisplay: Set[String],
                     labelSize: Double = 0.6,
                     ylab: String =
                       "depletion <-- Essential in cell lines --> enrichment",
                     doLegend: Boolean = true) = {
      val filteredFrame1 = frame.row(colored: _*)
      val numTest = math.log10(0.05 / (frame.numRows * 4)) * -1

      val colors: Map[String, Double] =
        legend
          .map(x => keggH.get(x).getOrElse("Unknown"))
          .toSeq
          .distinct
          .sorted
          .zipWithIndex
          .map(x => x._1 -> (x._2 + 1d))
          .toMap

      def caseTransform(s: String) =
        s.split("_")
          .map(x => x.head.toUpper +: x.tail.toLowerCase)
          .mkString(" ")

      val filteredLayer1 =
        (filteredFrame1
           .rconcat(
             Frame(
               "color" -> Series(filteredFrame1.rowIx.toSeq
                                   .map(x =>
                                     colors(keggH.get(x).getOrElse("Unknown")))
                                   .toVec,
                                 filteredFrame1.rowIx)))
           .mapRowIndex(rx =>
             if (labelsToDisplay.contains(rx))
               caseTransform(rx.stripPrefix("KEGG_"))
             else ""),
         point(size = 3d,
               labelText = true,
               labelFontSize = labelSize fts,
               color = DiscreteColors(colors.size),
               shapeCol = 4,
               sizeCol = 4))

      xyplot(
        frame -> point(size = 1.4d, color = Color.gray3),
        // filteredLayer2,
        filteredLayer1
      )(origin = Some(0d -> 0d),
        frame = false,
        xlab = "depletion <-- Essential in humans --> enrichment",
        ylab = ylab,
        xNumTicks = 0,
        yNumTicks = 0,
        xnames = Seq(numTest -> "", -1 * numTest -> ""),
        ynames = Seq(numTest -> "", -1 * numTest -> ""),
        xAxisMargin = 0d,
        yAxisMargin = 0d,
        xLabFontSize = .5 fts,
        yLabFontSize = .5 fts,
        xgrid = true,
        ygrid = true,
        // xlim = Some(-60d -> 60d),
        // ylim = Some(-60d -> 60d),
        legendFontSize = 0.35 fts,
        xCustomGrid = true,
        yCustomGrid = true,
        // xLineWidthFraction = 0.5,
        // yLineWidthFraction = 0.65,
        // yLineStartFraction = 0.25,
        // xLineStartFraction = 0.35,
        xTickLength = 0d fts,
        yTickLength = 0d fts,
        legendLayout = RelativeToFirst(10d, 10d),
        rightPadding = 35d,
        topPadding = 5d,
        extraLegend =
          if (doLegend)
            colors.toSeq.map(
              x =>
                x._1 -> PointLegend(shape = Shape.circle(1),
                                    color = DiscreteColors(colors.size)(x._2)))
          else Nil)

    }

    val keggLabelsInVitro = Set("KEGG_AMINOACYL_TRNA_BIOSYNTHESIS",
                                "KEGG_RIBOSOME",
                                "KEGG_DNA_REPLICATION",
                                "KEGG_PROTEASOME",
                                // "KEGG_CELL_CYCLE",
                                "KEGG_SPLICEOSOME",
                                // "KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
                                // "KEGG_MAPK_SIGNALING_PATHWAY",
                                // "KEGG_FOCAL_ADHESION",
                                // "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                                // "KEGG_OOCYTE_MEIOSIS",
                                "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
                                "KEGG_CALCIUM_SIGNALING_PATHWAY",
                                // "KEGG_PATHWAYS_IN_CANCER",
                                // "KEGG_LONG_TERM_POTENTIATION",
                                "KEGG_AXON_GUIDANCE",
                                "KEGG_RNA_DEGRADATION",
                                "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
                                "KEGG_BASAL_TRANSCRIPTION_FACTORS",
                                // "KEGG_MELANOGENESIS",
                                "KEGG_OLFACTORY_TRANSDUCTION")

    val keggLabelsMgi = Set( //"KEGG_AMINOACYL_TRNA_BIOSYNTHESIS",
                            "KEGG_RIBOSOME",
                            "KEGG_DNA_REPLICATION",
                            "KEGG_PROTEASOME",
                            // "KEGG_CELL_CYCLE",
                            "KEGG_SPLICEOSOME",
                            "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",
                            // "KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
                            // "KEGG_MAPK_SIGNALING_PATHWAY",
                            "KEGG_FOCAL_ADHESION",
                            // "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                            // "KEGG_OOCYTE_MEIOSIS",
                            // "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
                            // "KEGG_CALCIUM_SIGNALING_PATHWAY",
                            "KEGG_PATHWAYS_IN_CANCER",
                            // "KEGG_LONG_TERM_POTENTIATION",
                            "KEGG_AXON_GUIDANCE",
                            // "KEGG_RNA_DEGRADATION",
                            // "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
                            // "KEGG_BASAL_TRANSCRIPTION_FACTORS",
                            // "KEGG_MELANOGENESIS",
                            "KEGG_OLFACTORY_TRANSDUCTION")

    val gobppathwayplot = {
      val (frame, colors) =
        pathwayplot1(invivoEssentialGenes,
                     invitroEssentialGenes,
                     goBP,
                     table.rowIx.toSeq.toSet)
      pathwayplot2(frame, colors, colors, Set())
    }
    pdfToFile(new java.io.File("overrep_gobp.pdf"), gobppathwayplot, 1000)
    pngToFile(new java.io.File("overrep_gobp.png"), gobppathwayplot, 2000)

    val (invitropathwayplot, mgipathwayplot) = {
      val (frame1, colors1) =
        pathwayplot1(invivoEssentialGenes,
                     invitroEssentialGenes,
                     kegg,
                     table.rowIx.toSeq.toSet)
      val (frame2, colors2) =
        pathwayplot1(
          invivoEssentialGenes,
          dickinsonEssential,
          kegg,
          table.rowIx.toSeq.toSet & (dickinsonEssential ++ dickinsonNonEssential))
      val invitro =
        pathwayplot2(frame1,
                     colors1,
                     (colors1 ++ colors2).distinct,
                     keggLabelsMgi)
      val mgi = pathwayplot2(
        frame2,
        colors2,
        (colors1 ++ colors2).distinct,
        keggLabelsMgi,
        ylab = "depletion <-- Essential in mice --> enrichment",
        doLegend = false)
      (invitro, mgi)
    }
    val compositePathwayPlot =
      group(
        group(invitropathwayplot,
              AlignTo.topLeftCorner(TextBox("B"), invitropathwayplot.bounds),
              FreeLayout),
        group(mgipathwayplot,
              AlignTo.topLeftCorner(TextBox("C"), mgipathwayplot.bounds),
              FreeLayout),
        TableLayout(2))
    pdfToFile(new java.io.File("overrep_composite.pdf"),
              compositePathwayPlot,
              1000)
    pngToFile(new java.io.File("overrep_composite.png"),
              compositePathwayPlot,
              2000)

    pdfToFile(new java.io.File("overrep_kegg.pdf"), invitropathwayplot, 1000)
    pngToFile(new java.io.File("overrep_kegg.png"), invitropathwayplot, 2000)

    pdfToFile(new java.io.File("overrep_kegg_mgi.pdf"), mgipathwayplot, 1000)
    pngToFile(new java.io.File("overrep_kegg_mgi.png"), mgipathwayplot, 2000)

    val invivoButNotInVitro = invivoEssentialGenes & invitroNotEssentialGenes
    val invitroButNotInVivo = invitroEssentialGenes & invivoNotEssentialGenes
    val invivoAndInvitro = invivoEssentialGenes & invitroEssentialGenes
    val neigherInvivoAndInvitro = invivoNotEssentialGenes & invitroNotEssentialGenes

    println(" CELL LINES ")
    println("in vivo: " + invivoEssentialGenes.size)
    println("in vitro: " + invitroEssentialGenes.size)
    println("in vivo & in vitro : " + invivoAndInvitro.size)
    println(
      "in vivo &~ in vitro : " + (invivoEssentialGenes &~ invitroEssentialGenes).size)
    println(
      "in vitro &~ in vivo : " + (invitroEssentialGenes &~ invivoEssentialGenes).size)

    println(
      "in vitro &~ in vivo high power : " + ((invitroEssentialGenes &~ invivoEssentialGenes) & phiPower
        .filter(_ > 0.2)
        .index
        .toSeq
        .toSet).size)

    val invitroOnlyHighPower = (invitroEssentialGenes &~ invivoEssentialGenes) & phiPower
        .filter(_ > 0.2)
        .index
        .toSeq
        .toSet

    println(
      " in vitro &~ invivoEssentialGenes high power enrichment in recessive diseases: " + (overrepresentation
        .enrichmentTest(
          total = table.rowIx.toSeq.toSet.size,
          marked = (recessiveDiseaseGenes & (table.rowIx.toSeq.toSet)).size,
          draws = (invitroOnlyHighPower).size,
          markedDraws = (invitroOnlyHighPower & recessiveDiseaseGenes).size)))

    val highpower = phiPower.filter(_ > 0.2).index.toSeq.toSet

    println(" CELL LINE + MGI + IN VIVO")
    println(
      "in vivo + in vitro + mgi: " + (invivoEssentialGenes ++ invitroEssentialGenes ++ dickinsonEssential).size)
    println("in vivo: " + invivoEssentialGenes.size)
    println("in vitro: " + invitroEssentialGenes.size)
    println("mgi " + dickinsonEssential.size)
    println(
      "mgi & in vivo & in vitro: " + (dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes).size)
    println(
      "mgi &~ in vivo &~ in vitro: " + ((dickinsonEssential &~ invivoEssentialGenes) &~ invitroEssentialGenes).size)
    println(
      "in vivo &~ mgi &~ in vitro: " + ((invivoEssentialGenes &~ dickinsonEssential) &~ invitroEssentialGenes).size)
    println(
      "in vitro &~ mgi &~ in vivo: " + ((invitroEssentialGenes &~ dickinsonEssential) &~ invivoEssentialGenes).size)
    println(
      "(in vitro & mgi) &~ in vivo: " + ((invitroEssentialGenes & dickinsonEssential) &~ invivoEssentialGenes).size)
    println(
      "(in vivo & mgi) &~ in vitro: " + ((invivoEssentialGenes & dickinsonEssential) &~ invitroEssentialGenes).size)
    println(
      "(in vivo & in vitro) &~ mgi: " + ((invivoEssentialGenes & invitroEssentialGenes) &~ dickinsonEssential).size)

    println(
      "mgi &~ in vitro : " + (dickinsonEssential &~ invitroEssentialGenes).size)

    println(
      "in vivo &~ in vitro : " + (invivoEssentialGenes &~ invitroEssentialGenes).size)
    println(
      "in vitro &~ in vivo : " + (invitroEssentialGenes &~ invivoEssentialGenes).size)

    println(
      "in vitro & 0lof: " + (invitroEssentialGenes & lofcount
        .filter(_ == 0)
        .index
        .toSeq
        .toSet).size)
    println(
      "mgi & 0lof: " + (dickinsonEssential & lofcount
        .filter(_ == 0)
        .index
        .toSeq
        .toSet).size)

    println(
      "mgi & in vivo & in vitro: " + (dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes)
        .map(ensg2symbol))

    println(
      "mgi & in vivo & in vitro: \n" + (dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes)
        .map(x => ensg2symbol(x) -> unknowngenes(x))
        .toSeq
        .sortBy(_._1)
        .map(x => x._1 + "\t" + x._2)
        .mkString("\n"))

    {
      println(
        "CORE essential genes, network, dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes")
      val ((enrichment: Seq[(String, Double, Set[String])]), _) = overrep(
        (dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes),
        goBPShort,
        (dickinsonEssential ++ dickinsonNonEssential) & table.rowIx.toSeq.toSet)

      println(
        "Core essential genes: " + (dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes).size)

      println(
        "Core essential genes: \n" + (dickinsonEssential & invivoEssentialGenes & invitroEssentialGenes)
          .mkString("\n"))

      val manualSuperGroups = List(
        List("SPLICING",
             "SPLICEOSOME_ASSEMBLY",
             "MRNA_SPLICE_SITE_SELECTION",
             "RNA_SPLICINGVIA_TRANSESTERIFICATION_REACTIONS",
             "RNA_SPLICING",
             "RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS_AND_ASSEMBLY",
             "RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS_AND_ASSEMBLY",
             "PROTEIN_RNA_COMPLEX_ASSEMBLY",
             "MRNA_PROCESSING_GO_0006397",
             "MRNA_METABOLIC_PROCESS"),
        List("NUCLEAR_TRANSPORT",
             "NUCLEOBASENUCLEOSIDENUCLEOTIDE_AND_NUCLEIC_ACID_TRANSPORT",
             "NUCLEAR_TRANSPORT",
             "NUCLEAR_EXPORT",
             "RNA_EXPORT_FROM_NUCLEUS",
             "RNA_CATABOLIC_PROCESS",
             "NUCLEOCYTOPLASMIC_TRANSPORT"),
        List(
          "TRANSCRIPTION/EPIGENETICS/NUCLEOSOME",
          "NUCLEOSOME_ASSEMBLY",
          "TRANSCRIPTION_INITIATION",
          "PROTEIN_DNA_COMPLEX_ASSEMBLY",
          "CHROMATIN_REMODELING",
          "CHROMATIN_MODIFICATION",
          "ESTABLISHMENT_AND_OR_MAINTENANCE_OF_CHROMATIN_ARCHITECTURE",
          "ESTABLISHMENT_AND_OR_MAINTENANCE_OF_CHROMATIN_ARCHITECTURE",
          "REGULATION_OF_GENE_SPECIFIC_TRANSCRIPTION",
          "NEGATIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER",
          "REGULATION_OF_GENE_EXPRESSION_EPIGENETIC",
          "GENE_SILENCING"),
        List("MICROTUBULE/MITOSIS",
             "MICROTUBULE_BASED_PROCESS",
             "MICROTUBULE_CYTOSKELETON_ORGANIZATION_AND_BIOGENESIS",
             "MITOTIC_SPINDLE_ORGANIZATION_AND_BIOGENESIS",
             "MITOSIS",
             "M_PHASE_OF_MITOTIC_CELL_CYCLE"),
        List("CELL_CYCLE",
             "NEGATIVE_REGULATION_OF_CELL_CYCLE",
             "CELL_CYCLE_ARREST_GO_0007050",
             "MEIOTIC_CELL_CYCLE",
             "DOUBLE_STRAND_BREAK_REPAIR",
             "DNA_RECOMBINATION",
             "INDUCTION_OF_APOPTOSIS_BY_INTRACELLULAR_SIGNALS"))
        .flatMap(l => l.tail.map(x => x -> l.head))
        .toMap

      val bhThreshold =
        highestSignificantPValueByFDR(0.1, enrichment.map(_._2))
      println("bh threshold: " + bhThreshold)

      val edgeListFromString =
        readStringTsv(stringCoreTsv, symbol2ensg, hgncTable)

      println("Edges: " + edgeListFromString.size)

      println(
        "Nodes: " + edgeListFromString
          .flatMap(x => List(x._1, x._2))
          .distinct
          .size)

      val significantEnrichment =
        enrichment.filter(_._2 < (bhThreshold))

      val genesWithSuperGroups = significantEnrichment.flatMap {
        case (go, p, genes) =>
          genes.map(g => (g, go, manualSuperGroups.get(go).getOrElse("Other")))
      }.map(x => x._1 -> (x._2 -> x._3)).toMap

      println(significantEnrichment.mkString("\n"))

      writeToFile("annotated_core_string_edges.csv",
                  "g1,g2,v\n" + edgeListFromString.map {
                    case (g1, g2, v) =>
                      List(ensg2symbol(g1), ensg2symbol(g2), v).mkString(",")
                  }.mkString("\n"))

      writeToFile(
        "annotated_core_string_nodes.csv",
        "g,go,sg\n" + edgeListFromString
          .flatMap(x => List(x._1, x._2))
          .map { g =>
            ensg2symbol(g) + "," + genesWithSuperGroups
              .get(g)
              .map(_._1)
              .getOrElse("ns") + "," +
              genesWithSuperGroups.get(g).map(_._2).getOrElse("ns")
          }
          .mkString("\n")
      )

    }

    {
      println("Pairwise comparison with Odds-Ratio")
      val contingency = percentiles.colIx.toSeq.flatMap { cx1 =>
        percentiles.colIx.toSeq.map { cx2 =>
          val joined = percentiles.col(cx1, cx2).rdropNA > percentileCutoff

          val both =
            joined
              .rfilter(series => series.get(cx1).get && series.get(cx2).get)
              .numRows
          val neither =
            joined
              .rfilter(series => !series.get(cx1).get && !series.get(cx2).get)
              .numRows
          val a1 =
            joined
              .rfilter(series => series.get(cx1).get && !series.get(cx2).get)
              .numRows
          val a2 =
            joined
              .rfilter(series => !series.get(cx1).get && series.get(cx2).get)
              .numRows

          (cx1, cx2, ContingencyTable2x2(both, a2, a1, neither))

        }
      }

      val maxC = contingency.map(x => or(x._3)).filterNot(_.isInfinite).max

      def replaceInf1(x: Double) = if (x.isInfinite) maxC else x

      val orplot = Heatmap
        .fromColumns(
          frame = contingency
            .map(x => (x._1, x._2, math.log10(replaceInf1(or(x._3)))))
            .toFrame,
          reorderRows = true,
          euclidean = true,
          valueText = false,
          // zlim = Some(-1d -> 1d),
          colormap = RedBlue(-1, 1, 0)
        )
        ._1
      pdfToFile(new java.io.File("or.pdf"), orplot, 1000)
      pngToFile(new java.io.File("or.png"), orplot, 2000)

    }

    (table, percentiles)

  }
}
