/**
  * Represents a 2x2 Contingency table
  * @param a11 First row first column
  * @param a12 First row second column
  * @param a21 Second row first column
  * @param a22 Second row second column
  */
// ContingencyTable2x2(both, a2, a1, neither)
case class ContingencyTable2x2(a11: Int, a12: Int, a21: Int, a22: Int) {
  def sum = a11 + a12 + a21 + a22
  def rowMargin1 = a11 + a12
  def rowMargin2 = a21 + a22
  def columnMargin1 = a11 + a21
  def columnMargin2 = a12 + a22
}

object Main extends App {
  Runner.run

}
