import overrepresentation._

object OverRep {

  /**
    * FDR Benjamini-Cochberg procedure
    *
    * ref: http://udel.edu/~mcdonald/statmultcomp.html
    */
  def highestSignificantPValueByFDR(Q: Double,
                                    pValues: Iterable[Double]): Double = {
    val Qpm = Q / pValues.size.toDouble
    pValues.toSeq.sorted.zipWithIndex.reverse
      .find(x => x._1 <= (x._2 + 1).toDouble * Qpm)
      .map(_._1) match {
      case Some(x) => x
      case None => 0.0
    }
  }

  def overrep(set: Set[String],
              pathways: Seq[(String, Set[String])],
              bg: Set[String]) = {

    (pathways
       .map(_._1)
       .zip(enrichmentTestWithBackground(pathways.map(_._2), set, bg))
       .map(x => (x._1, x._2._1, set & pathways.find(_._1 == x._1).get._2)),
     pathways
       .map(_._1)
       .zip(depletionTestWithBackground(pathways.map(_._2), set, bg))
       .map(x => (x._1, x._2._1, set & pathways.find(_._1 == x._1).get._2)))
  }

  def enrichmentTestWithBackground[T](aprioriSets: Seq[Set[T]],
                                      targetSet: Set[T],
                                      bg: Set[T]): Seq[(Double, Counts)] = {
    val bgs = bg.size
    val tgbg = (bg & targetSet)
    val tgbgs = tgbg.size
    aprioriSets.map { ap =>
      val md = (tgbg & ap).size
      val c =
        Counts(total = bgs,
               marked = (ap & bg).size,
               draws = tgbgs,
               markedDraws = md)
      val p = overrepresentation.enrichmentTest(c)
      (p, c)
    }
  }

  def depletionTestWithBackground[T](aprioriSets: Seq[Set[T]],
                                     targetSet: Set[T],
                                     bg: Set[T]): Seq[(Double, Counts)] = {

    val bgs = bg.size
    val tgbg = (bg & targetSet)
    val tgbgs = tgbg.size
    aprioriSets.map { ap =>
      val md = (tgbg & ap).size
      val c =
        Counts(total = bgs,
               marked = (ap & bg).size,
               draws = tgbgs,
               markedDraws = md)
      val p = overrepresentation.depletionTest(c)
      (p, c)
    }
  }

}
