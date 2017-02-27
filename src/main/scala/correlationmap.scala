import org.saddle._
import org.saddle.linalg._
import org.nspl._
import org.nspl.saddle._

object CorrelationPlot {

  def fromColumns[RX <: AnyRef: ST: ORD, CX <: AnyRef: ST: ORD](
      f: Frame[RX, CX, Double],
      main: String = "",
      xlab: String = "",
      ylab: String = "",
      xLabFontSize: Option[RelFontSize] = None,
      yLabFontSize: Option[RelFontSize] = None,
      mainFontSize: RelFontSize = 1 fts,
      colormap: Colormap = RedBlue(-1, 1, 0)) = {

    val vecs = f.toColSeq.map {
      case (cx, series) =>
        (cx, series.toVec)
    }

    val ar = Array.ofDim[Double](vecs.size * vecs.size)

    var i = 0
    var j = 0
    while (i < vecs.size) {
      while (j <= i) {
        val (c1, v1) = vecs(i)
        val (c2, v2) = vecs(j)

        val missing = (v1.toSeq.zipWithIndex.filter(_._1.isNaN).map(_._2) ++
          v2.toSeq.zipWithIndex.filter(_._1.isNaN).map(_._2)).toArray

        val v1a = v1.without(missing)
        val v2a = v2.without(missing)

        val v1ad = v1a.demeaned
        val v2ad = v2a.demeaned
        val s1 = v1a.stdev
        val s2 = v2a.stdev

        val cov = v1ad vv v2ad * (1d / (v1ad.length - 1))
        val r = cov / (s1 * s2)
        ar(i * vecs.size + j) = r
        ar(j * vecs.size + i) = r
        j += 1
      }
      j = 0
      i += 1
    }

    val f2 = Frame(Mat(vecs.size, vecs.size, ar),
                   Index(vecs.map(_._1): _*),
                   Index(vecs.map(_._1): _*))

    val plot = Heatmap
      .fromColumns(
        frame = f2,
        reorderRows = true,
        euclidean = false,
        main = main,
        xlab = xlab,
        ylab = ylab,
        xLabFontSize = xLabFontSize,
        yLabFontSize = yLabFontSize,
        mainFontSize = mainFontSize,
        colormap = colormap,
        valueText = false,
        zlim = Some(-1d -> 1d)
      )
      ._1

    (plot, f2)

  }

}
