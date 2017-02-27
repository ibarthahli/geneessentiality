package org.nspl

import data._

case class DataElem(
    data: DataSource,
    xAxis: Axis,
    yAxis: Axis,
    renderers: Seq[DataRenderer],
    bounds: Bounds,
    tx: AffineTransform = AffineTransform.identity
) extends Renderable[DataElem] {
  def transform(tx: Bounds => AffineTransform) = this.copy(tx = tx(bounds).concat(this.tx))
}

trait Plots {

  type XYPlotArea = Elems5[ElemList[ShapeElem], ElemList[ShapeElem], ElemList[DataElem], Elems2[AxisElem, AxisElem], ShapeElem]

  type Figure = Elems3[ShapeElem, Elems3[ShapeElem, Elems2[Elems2[Elems2[XYPlotArea, TextBox], TextBox], TextBox], ShapeElem], ShapeElem]

  def xyplotarea(
    data: Seq[(DataSource, List[DataRenderer])],
    xAxisSetting: AxisSettings = AxisSettings(LinearAxisFactory),
    yAxisSetting: AxisSettings = AxisSettings(LinearAxisFactory),
    origin: Option[Point] = None,
    xlim: Option[(Double, Double)] = None,
    ylim: Option[(Double, Double)] = None,
    xAxisMargin: Double = 0.05,
    yAxisMargin: Double = 0.05,
    xCol: Int = 0,
    yCol: Int = 1,
    xgrid: Boolean = true,
    ygrid: Boolean = true,
    frame: Boolean = true,
    boundsData: Seq[DataSource] = Nil,
    xCustomGrid: Boolean = false,
    yCustomGrid: Boolean = false,
    main: String = "",
    mainFontSize: RelFontSize = 1.2 fts,
    mainDistance: RelFontSize = 1.2 fts,
    xlab: String = "",
    xlabFontSize: RelFontSize = 1.0 fts,
    xlabDistance: RelFontSize = 1.0 fts,
    xlabAlignment: Alignment = Center,
    ylab: String = "",
    ylabFontSize: RelFontSize = 1.0 fts,
    ylabDistance: RelFontSize = 1.0 fts,
    ylabAlignment: Alignment = Center,
    topPadding: Double = 2d,
    bottomPadding: Double = 0d,
    leftPadding: Double = 0d,
    rightPadding: Double = 2d
  ): Figure = {

    val xMinMax =
      if (boundsData.isEmpty) data.map(_._1.columnMinMax(xCol))
      else boundsData.map(_.columnMinMax(xCol))

    val yMinMax =
      if (boundsData.isEmpty) data.map(_._1.columnMinMax(yCol))
      else boundsData.map(_.columnMinMax(yCol))

    val xLimMin = xlim.map(_._1).filterNot(_.isNaN)
    val xLimMax = xlim.map(_._2).filterNot(_.isNaN)

    val yLimMin = ylim.map(_._1).filterNot(_.isNaN)
    val yLimMax = ylim.map(_._2).filterNot(_.isNaN)

    val dataXMin = if (xLimMin.isDefined) 0.0 else xMinMax.map(_.min).min
    val dataXMax = if (xLimMax.isDefined) 0.0 else xMinMax.map(_.max).max
    val dataYMin = if (yLimMin.isDefined) 0.0 else yMinMax.map(_.min).min
    val dataYMax = if (yLimMax.isDefined) 0.0 else yMinMax.map(_.max).max

    val xMin = math.min(xLimMin.getOrElse {
      dataXMin - xAxisMargin * (dataXMax - dataXMin)
    }, origin.map(_.x).getOrElse(Double.MaxValue))

    val xMax = xLimMax.getOrElse {
      dataXMax + xAxisMargin * (dataXMax - dataXMin)
    }

    val yMin = math.min(yLimMin.getOrElse {
      dataYMin - yAxisMargin * (dataYMax - dataYMin)
    }, origin.map(_.y).getOrElse(Double.MaxValue))

    val yMax = yLimMax.getOrElse {
      dataYMax + yAxisMargin * (dataYMax - dataYMin)
    }

    val xAxis = xAxisSetting.axisFactory.make(xMin, xMax, xAxisSetting.width, true)
    val yAxis = yAxisSetting.axisFactory.make(yMin, yMax, yAxisSetting.width, false)

    val xMinV = xAxis.worldToView(xMin)
    val xMaxV = xAxis.worldToView(xMax)
    val yMinV = yAxis.worldToView(yMin)
    val yMaxV = yAxis.worldToView(yMax)

    val yAxisViewMin = yAxisSetting.lineStartFraction * yAxis.width
    val yAxisViewMax = yAxisViewMin + yAxisSetting.width * yAxisSetting.lineLengthFraction

    val xAxisViewMin = xAxisSetting.lineStartFraction * xAxis.width
    val xAxisViewMax = xAxisViewMin + xAxisSetting.width * xAxisSetting.lineLengthFraction

    val originWX1 = origin.map { origin =>
      xlim.map {
        case (a, b) =>
          if (origin.x < a) a
          else if (origin.x > b) b
          else origin.x
      } getOrElse origin.x
    } getOrElse xMin

    val originWY1 = origin.map { origin =>
      ylim.map {
        case (a, b) =>
          if (origin.y < a) a
          else if (origin.y > b) b
          else origin.y
      } getOrElse origin.y
    } getOrElse yMin

    val noXTick = if (origin.isEmpty) Nil else List(originWX1)
    val noYTick = if (origin.isEmpty) Nil else List(originWY1)

    val (xMajorTicks, xCustomTicks, xAxisElem) = xAxisSetting.renderable(xAxis, noXTick)
    val (yMajorTicks, yCustomTicks, yAxisElem) = yAxisSetting.renderable(
      yAxis,
      noYTick
    )

    val originX = xAxis.worldToView(originWX1)
    val originY = yAxis.worldToView(originWY1)

    val axes = group(
      translate(xAxisElem, 0, originY),
      translate(yAxisElem, originX, 0),
      FreeLayout
    )

    val dataelem = sequence(data.toList.map {
      case (ds, drs) =>
        DataElem(ds, xAxis, yAxis, drs, axes.bounds, AffineTransform.identity)
    })

    val xgridPoints = if (xCustomGrid) (xMajorTicks ++ xCustomTicks).distinct else xMajorTicks

    val ygridPoints = if (yCustomGrid) (yMajorTicks ++ yCustomTicks).distinct else yMajorTicks

    val xgridElem = sequence(xgridPoints map { w =>
      val v = xAxis.worldToView(w)
      ShapeElem(
        Shape.line(Point(v, yAxisViewMin), Point(v, yAxisViewMax)),
        stroke = if (xgrid) Some(Stroke(1d)) else None,
        strokeColor = Color.gray5
      )
    })

    val frameStroke = if (frame) Some(Stroke(1d)) else None
    val frameElem =
      ShapeElem(
        Shape.rectangle(xMinV, yMaxV, xMaxV - xMinV, math.abs(yMinV - yMaxV)), stroke = frameStroke, fill = Color.transparent
      )

    val ygridElem = sequence(ygridPoints map { w =>
      val v = yAxis.worldToView(w)
      ShapeElem(
        Shape.line(
          Point(xAxisViewMin, v),
          Point(xAxisViewMax, v)
        ),
        stroke = if (ygrid) Some(Stroke(1d)) else None,
        strokeColor = Color.gray5
      )
    })

    val renderedPlot = group(xgridElem, ygridElem, dataelem, axes, frameElem, FreeLayout)

    val mainBox = TextBox(main, fontSize = mainFontSize, width = Some(frameElem.bounds.w))
    val xlabBox = TextBox(xlab, fontSize = xlabFontSize, width = Some(frameElem.bounds.w))
    val ylabBox = TextBox(ylab, fontSize = ylabFontSize, width = Some(frameElem.bounds.h))

    val padTop = ShapeElem(
      shape = Shape.circle(topPadding),
      fill = Color.transparent,
      strokeColor = Color.transparent,
      stroke = None
    )

    val padBottom = ShapeElem(
      shape = Shape.circle(bottomPadding),
      fill = Color.transparent,
      strokeColor = Color.transparent,
      stroke = None
    )

    val padLeft = ShapeElem(
      shape = Shape.circle(leftPadding),
      fill = Color.transparent,
      strokeColor = Color.transparent,
      stroke = None
    )

    val padRight = ShapeElem(
      shape = Shape.circle(rightPadding),
      fill = Color.transparent,
      strokeColor = Color.transparent,
      stroke = None
    )

    val withHorizontalLabels = zgroup(
      (zgroup(
        (renderedPlot, 1),
        (AlignTo.horizontalCenter(mainBox, frameElem.bounds), 0),
        VerticalStack(NoAlignment, gap = mainDistance)
      ), 0),
      (AlignTo.horizontal(xlabBox, frameElem.bounds, xlabAlignment), 1),
      VerticalStack(NoAlignment, xlabDistance)
    )

    val movedFrame = withHorizontalLabels.m1.m1.m5

    val plotWithAxisLabels =
      zgroup(
        (withHorizontalLabels, 1),
        (AlignTo.vertical(rotate(ylabBox, 0.5 * math.Pi), movedFrame.bounds, ylabAlignment), 0),

        HorizontalStack(NoAlignment, ylabDistance)
      )

    group(padTop, group(padLeft, plotWithAxisLabels, padRight, HorizontalStack(Center, 0d)), padBottom, VerticalStack(Center, 0d))

  }

  sealed trait LegendElem
  case class PointLegend(shape: Shape, color: Color) extends LegendElem
  case class LineLegend(stroke: Stroke, color: Color) extends LegendElem

  type Legend = ElemList[Elems2[ShapeElem, TextBox]]

  def legend(
    entries: List[(String, LegendElem)],
    fontSize: RelFontSize = 1.0 fts,
    width: RelFontSize = 30 fts
  ): Legend = {
    sequence(entries.map {
      case (text, elem) =>
        val elem1 = elem match {
          case PointLegend(s, c) => fitToBounds(ShapeElem(s, fill = c), Bounds(0, 0, fontSize, fontSize))
          case LineLegend(s, c) =>
            fitToBounds(
              ShapeElem(
                Shape.line(Point(0, 0), Point(fontSize, 0)),
                strokeColor = c,
                stroke = Some(s)
              ),
              Bounds(0, 0, fontSize, fontSize)
            )
        }
        group(
          elem1,
          TextBox(text, fontSize = fontSize, width = Some(width)),
          HorizontalStack(Center, fontSize)
        )
    }, VerticalStack(Left))
  }

  type HeatmapLegend = Elems2[ElemList[ShapeElem], Elems3[ShapeElem, ElemList[Elems2[ShapeElem, TextBox]], ElemList[ShapeElem]]]

  def heatmapLegend(
    min: Double,
    max: Double,
    color: Colormap = HeatMapColors(0d, 1d),
    fontSize: RelFontSize = 1.0 fts,
    width: RelFontSize = 10 fts,
    height: RelFontSize = 1 fts
  ): HeatmapLegend = {

    val color1 = color.withRange(min, max)

    val axisSettings = AxisSettings(
      LinearAxisFactory,
      fontSize = fontSize,
      width = width,
      numTicks = 2,
      tickAlignment = 1
    )
    val axis = axisSettings.axisFactory.make(min, max, width, false)

    val (majorTicks, _, axisElem) =
      axisSettings.renderable(
        axis
      )

    val n = 500
    val space = (max - min) / n.toDouble
    val ticks = sequence(
      ((0 until n toList) map { i =>
        val world = axis.min + i * space
        val view = axis.worldToView(world)
        ShapeElem(
          Shape.line(Point(1d, view), Point(height, view)),
          stroke = Some(Stroke(2d)),
          strokeColor = color1(world)
        )
      })
    )

    group(ticks, axisElem, FreeLayout)

  }

}
