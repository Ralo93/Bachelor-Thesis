====== ALEVIN - ALgorithms for Embedding Virtual Networks ======

Recommended IDE:
Eclipse version 3.5 or higher

===== Requirements =====

  * JUNG 2.0.1 or newer, http://jung.sf.net
  * JUnit 4.5 or newer, http://www.junit.org
  * GLPK for Java 1.0.x or newer, http://glpk-java.sourceforge.net/
  * Apache commons CLI 1.2 or newer, https://commons.apache.org/proper/commons-cli/
  * Sun/Oracle JRE 1.6.x or newer

Note that some of these dependencies may also be available via your Operating
System (i.e., as a system wide installation). In that case, modify the
respective entries in ".classpath" to reflect the system wide location of the
respective libraries. Alternatively, in *nix systems it should be possible to
create symlinks in the lib/ folder to point to the right JAR file.

===== Installation =====

  * Run the ANT script build.xml (either from console or from Eclipse)
  * This will create the following files in the base directory
    * ALEVIN.jar
    * ALEVIN-src.jar
    * JUNG2.jar

!! Actually, automatic download of libraries is currently broken! Use the URLs
!! provided above to download the JAR files manually and put them in the lib/
!! folder.


===== Running ALEVIN =====

After running the build script, there are two options:
  * Inside of Eclipse, you can now simply run the vnreal.Main class
  * To run ALEVIN from JAR files, copy the JAR files created by the
    installation script to your destination folder, change into this directory
    and run
      java -jar ALEVIN.jar  
  * For running a large number of tests without GUI have a look at the files in
    src/tests/generatorTests/

===== Features =====

  * An XML-based scenario import and export, including
    * Topology data
    * Resources (substrate)
    * Demands (virtual networks)
  * Graphical user interface (GUI) based on MuLaViTo providing
    * QuickSearchBar
    * Display of resources and demands
    * Advanced graph editing and creation
      * Adding of resource/demands to network entities

===== Basic Structure =====

The Java package structure for this project looks like:
  * ''AUTHORS'' -- A list of contributors to this software and their affiliations
  * ''COPYING'' -- The license information (GPL, LGPL)
  * ''README'' -- This file
  * ''build.xml'' -- The ANT build script
  * ''src''
    * ''vnreal''
      * ''mulavito'' --- contains the mulavito source
      * ''algorithms'' -- All kinds of algorithms
      * ''constraints'' -- Properties belonging to both, requests as well as resources
      * ''demands'' -- All classes representing virtual demands
      * ''evaluations'' -- Evaluation metrics and utils
      * ''gui'' -- GUI elements to come
      * ''hiddenhopmapping'' -- Hidden hop mapping functions
      * ''io'' -- Import/export of topologies
      * ''mapping'' -- The data structure storing the mapping of requests to resources and vice versa
      * ''network'' -- Basic network elements (networks, nodes, links)
      * ''resources'' -- All classes representing resources in the substrate
    * ''tests'' -- JUnit tests (usually not exported with ordinary source code)
  * ''ILP-LP-Models'' -- (Integer) Linear Programs for GLPK with gplk-java, see User guide
  * ''tools'' -- Some scripts to evaluate results

  
===== Known Issues / Limitations =====
 
  * Each demand only occupies a single resource.
  * jsc.jar (http://www.jsc.nildram.co.uk/downloads/download.html) is used in
    some evaluation tools - due to licensing issues this cannot be included in
    the SF SVN branch. If necessary, please download manually
