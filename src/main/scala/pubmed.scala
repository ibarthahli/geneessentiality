import play.api.libs.ws.ahc._
import play.api.libs.ws._
import akka.actor._
import akka.stream.ActorMaterializer
import scala.concurrent.duration._
import scala.concurrent._

import scala.concurrent.ExecutionContext.Implicits.global

object Pubmed {
  def queryCount(term: String)(implicit ws: WSClient): Future[Int] =
    ws.url("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi")
      .withQueryString("db" -> "pubmed",
                       "retmode" -> "json",
                       "retmax" -> "1",
                       "term" -> term)
      .get
      .map { response =>
        try { (response.json \ "esearchresult" \ "count").as[String] } catch {
          case e: Exception => {
            println(response.body)
            throw e
          }
        }
      }
      .map(_.toInt)
      .recover {
        case e => {
          println(e)
          0
        }
      }

  def queryCountSync(terms: Seq[String]) = {
    implicit val as = ActorSystem()
    implicit val materializer = ActorMaterializer()
    implicit val ws = AhcWSClient(
      ahc.AhcWSClientConfigFactory.forClientConfig(
        WSClientConfig(
          connectionTimeout = 2.minutes,
          idleTimeout = 2.minutes,
          requestTimeout = 2.minutes,
          followRedirects = true,
          useProxyProperties = true,
          userAgent = None,
          compressionEnabled = false,
          ssl = ssl.SSLConfig(
            loose = ssl.SSLLooseConfig(acceptAnyCertificate = true))
        )
      ))

    val r = terms.map { term =>
      val r = Await.result(queryCount(term)(ws), atMost = Duration.Inf)
      println((term, r))
      (term, r)
    }

    ws.close
    as.shutdown
    r
  }
}
