scalaVersion := "2.11.8"

libraryDependencies ++= Seq(
  "io.github.pityka" %% "fileutils" % "1.0.0",
  "io.github.pityka" %% "stringsplit" % "1.0.0",
  "io.github.pityka" %% "overrepresentation" % "1.0.0",
  "io.github.pityka" %% "saddle-linalg" % "0.0.15",
  "org.scala-saddle" %% "saddle-core" % "1.3.4",
  "io.github.pityka" % "hierarchical-clustering-fork" % "1.0-5",
  "org.apache.commons" % "commons-math3" % "3.6.1",
  "com.typesafe.play" %% "play-ws" % "2.5.9"
).map(
  x =>
    x.exclude("log4j", "log4j")
      .exclude("commons-logging", "commons-logging")
      .excludeAll(ExclusionRule(organization = "ch.qos.logback"))
)

lazy val commonSettings = Seq(scalaVersion := "2.11.8")

lazy val core = project
  .in(file("nspl/core"))
  .settings(commonSettings)
  .settings(
    name := "nspl-core"
  )
  .enablePlugins(spray.boilerplate.BoilerplatePlugin)

lazy val awt = project
  .in(file("nspl/awt"))
  .settings(commonSettings)
  .settings(
    name := "nspl-awt",
    libraryDependencies += "de.erichseifert.vectorgraphics2d" % "VectorGraphics2D" % "0.11"
  )
  .dependsOn(core)

lazy val saddle = (project in file("nspl/saddle"))
  .settings(commonSettings)
  .settings(
    name := "nspl-saddle",
    libraryDependencies ++= Seq(
      "org.scala-saddle" %% "saddle-core" % "1.3.4" exclude ("com.googlecode.efficient-java-matrix-library", "ejml"),
      "com.googlecode.efficient-java-matrix-library" % "ejml" % "0.19" % "test",
      "org.scalatest" %% "scalatest" % "2.1.5" % "test"
    )
  )
  .dependsOn(core, awt)

lazy val root = project.in(file(".")).dependsOn(saddle, awt)

reformatOnCompileSettings
