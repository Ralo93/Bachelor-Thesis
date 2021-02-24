package vnreal.algorithms;

import java.sql.Time;
import org.jfree.chart.ChartFactory;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.Vector;
import java.util.stream.DoubleStream;

import javax.swing.JFrame;

import io.jenetics.internal.math.random;
import jsc.util.BitVector;
import mulavito.algorithms.AbstractAlgorithmStatus;
import tests.VNData;
import vnreal.algorithms.linkmapping.kShortestPathLinkMapping;
import vnreal.algorithms.linkmapping.kShortestPathLinkMappingTest;
import vnreal.algorithms.utils.Chromosome;
import vnreal.algorithms.utils.MiscelFunctions;
import vnreal.constraints.demands.AbstractDemand;
import vnreal.constraints.demands.BandwidthDemand;
import vnreal.constraints.demands.CpuDemand;
import vnreal.constraints.resources.AbstractResource;
import vnreal.constraints.resources.BandwidthResource;
import vnreal.constraints.resources.CpuResource;
import vnreal.evaluations.metrics.Cost;
import vnreal.gui.mapping.MappingTreeSelectionListener;
import vnreal.mapping.Mapping;
import vnreal.network.Network;
import vnreal.network.NetworkStack;
import vnreal.network.Node;
import vnreal.network.substrate.SubstrateLink;
import vnreal.network.substrate.SubstrateNetwork;
import vnreal.network.substrate.SubstrateNode;
import vnreal.network.virtual.VirtualLink;
import vnreal.network.virtual.VirtualNetwork;
import vnreal.network.virtual.VirtualNode;

import java.awt.BorderLayout;
import java.awt.color.*;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JDialog;
import javax.swing.JPanel;

public class GeneticAlgo extends AbstractSequentialAlgorithm<VirtualNetwork> {

	long starttime = 0;
	long endtime = 0;
	long runtime = 0;
	private final NetworkStack ns;
	private Iterator<? extends Network<?, ?, ?>> curNetIt;
	private Iterator<VirtualNetwork> curIt;
	private int processedVNetworks = 0;
	private int hostedNetworks = 0;
	// Koennen auch noch mehrere Metriken g
	private FileWriter writer = null;
	private int reinit = 0;
	private int broken = 0;
	private long timeLocal = 0;
	private int localSuccess = 0;
	private static int mutCheck = 0;

	public GeneticAlgo(NetworkStack ns) {
		this.ns = ns;
		this.curNetIt = ns.iterator();
//		this.curIt = null;
		// ns.getVirtuals().iterator();
	}

	public int getHosted() {

		return this.hostedNetworks;
	}

	@Override
	protected boolean preRun() {
		//
		try {
			writer = new FileWriter("10BasicBigTest10.csv");
		} catch (IOException e) {
			System.out.println("Bad writer!");
			e.printStackTrace();
		}
		this.starttime = System.currentTimeMillis();
		System.out.println("This is the GeneticAlgo Prerun at " + this.starttime);
		return true;

	}

	@Override
	protected void postRun() {
		try {
			writer.flush();
			writer.close();
		} catch (IOException e) {
			System.out.println("Closing failed in post run!");
			e.printStackTrace();
		}
		this.endtime = System.currentTimeMillis();
		this.runtime = this.endtime - this.starttime;
		double realruntime = this.runtime;
		System.out.println("This GeneticAlgo took " + realruntime + " Millisec");
		System.out.println("Under 110: " + this.broken);
		System.out.println("Local Time: " + String.valueOf(this.timeLocal));
		System.out.println("Local Sucesses: " + this.localSuccess);
		System.out.println("MutCheck: " + this.mutCheck);
	}

//	protected Chromosome easyInit(Chromosome c, SubstrateNetwork sn, VirtualNetwork vn, Map<VirtualNode, Integer> mapp) {
//		double[] hostvector = c.getHVector();

//		List<SubstrateNode> untouchedsubs = new LinkedList<SubstrateNode>();
//		for (SubstrateNode n : sn.getVertices()) {
//			n.setTouched(false);
//		}

	// for (int i = 0; i < hostvector.length; i++) {
	// if (((long)hostvector[i] == (long) n.getId())) {
	// n.setTouched(true);
//					System.out.println("SNode set to touched");
	// continue;
	// } else {
	// untouchedsubs.add(n);
	// continue;
	// }
	// }
//			System.out.println("A substrate node with touched: " + n.getTouched() + " and ID: " + n.getId());	
	// }

//		Queue<VirtualNode> pq = new LinkedList<VirtualNode>();
//		int i = 0;
//		for (VirtualNode vno : vn.getVertices()) {	
	// Nimm virtuellen knoten mit diesem index und pack ihn in PQ
//				pq.add(vno);
//				System.out.println("Virtual Node with ID: " + vno.getId() + " added to Queue!");
	// i++;
//		}

//	}

	// Algorithm 1 von Mi et al. 2012
	// stellt eine verbesserte Initialisierung dar, der parameter mapp sind
	// die indizes der VNodes des VNR
	// FIXME
	protected Chromosome improvedInit(Chromosome c, SubstrateNetwork sn, VirtualNetwork vn,
			Map<VirtualNode, Integer> mapp) {
		double[] hostvector = c.getHVector();
		boolean[] statevector = c.getSVector();
		// System.out.println("The VirtualNetwork: " + vn);

		// System.out.println("THIS IS ALGO1!");
		// for (double d : hostvector) {
		// System.out.println("The value of the hostvector for every entry: " +
		// String.valueOf(d));
		// }
		// for (boolean d : statevector) {
		// System.out.println("The value of the statevector for every entry: " + d);
		// }

		List<SubstrateNode> untouchedsubs = new LinkedList<SubstrateNode>();
		for (SubstrateNode n : sn.getVertices()) {

			n.setTouched(false);
			for (int i = 0; i < hostvector.length; i++) {
				if (((long) hostvector[i] == (long) n.getId()) && statevector[i]) {
					n.setTouched(true);
//					System.out.println("SNode set to touched");
					continue;
				} else {
//					System.out.println("SNode stays untouched");
					untouchedsubs.add(n);
					continue;
				}
			}
//			System.out.println("A substrate node with touched: " + n.getTouched() + " and ID: " + n.getId());	
		}

		// subnodes = (List<SubstrateNode>) sn.getVertices();
		Queue<VirtualNode> pq = new LinkedList<VirtualNode>();
		int i = 0;
		for (VirtualNode vno : vn.getVertices()) {
			System.out.println();
			if (!statevector[i]) {
				// Nimm virtuellen knoten mit diesem index und pack ihn in PQ
				pq.add(vno);
//				System.out.println("Virtual Node with ID: " + vno.getId() + " added to Queue!");
				i++;
				continue;
			} else {
				i++;
				continue;
			}
		}

//		System.out.println("PQ SIZE: " + pq.size() );
		Map<VirtualNode, Double> pqRanking;
		pqRanking = new LinkedHashMap<VirtualNode, Double>();
		Map<VirtualNode, Integer> originalPosition;
		originalPosition = new LinkedHashMap<VirtualNode, Integer>();

		Map<SubstrateNode, Double> untouchedRanking;
		untouchedRanking = new LinkedHashMap<SubstrateNode, Double>();
		for (SubstrateNode subnode : untouchedsubs) {
			untouchedRanking.put(subnode, MiscelFunctions.getAr(sn, subnode));
		}
		originalPositionLook(vn, pq, pqRanking, originalPosition, mapp);

		Map<VirtualNode, Double> sortedPQ = MiscelFunctions.sortByValue(pqRanking);
		// System.out.println("This is sorted now: " + sortedPQ);
		// Map<SubstrateNode, Double> sortedUntRank =
		// MiscelFunctions.sortByValue(untouchedRanking);
		// System.out.println("This is sorted now: " + sortedUntRank);
		// Siehe Schritt 4 in algo1
		Iterator<VirtualNode> nodesIt = sortedPQ.keySet().iterator();
//		Iterator<VirtualNode> origIt = originalPosition.keySet().iterator();
		// System.out.println("SUBSTRATE: " + sn);
		for (int in = 0; in < sortedPQ.size(); in++) {
			int originalPositionOfVNode = 0;
			// Map.Entry<VirtualNode, Double> entry = sortedPQ.entrySet().iterator().next();
			// make candidate list for every node in pq and map one of it to the hostvector
			List<SubstrateNode> candidates = new LinkedList<SubstrateNode>();
			VirtualNode vnode = nodesIt.next();
			originalPositionOfVNode = originalPosition.get(vnode);
			for (AbstractDemand dem : vnode) {

				if (dem instanceof CpuDemand) {

					candidates = findFulfillingNodes(dem);
//					System.out.println("Candidate Nodes for: " + vnode.getId() + " " + candidates);	
					double allprob = 0;
					double hostprob = 0;
					for (SubstrateNode subN : candidates) {

						double tmp = MiscelFunctions.getAr(sn, subN);

//						System.out.println("TMP: " + tmp);

						allprob = allprob + tmp;
//						System.out.println("Allprob nach der Miscelberechnung : " + allprob);
						// System.out.println("Allprob: " + allprob);
					}
					Double rand = Math.random();
//					System.out.println("Random: " + rand);
					Double[] wheel = new Double[candidates.size()];
//					System.out.println("Wheel Groesse: " + wheel.length);
					double counter = 0;
					for (int im = 0; im < candidates.size(); im++) {

						double percent = (MiscelFunctions.getAr(sn, candidates.get(im))) / allprob;
//						System.out.println("Miscel for the subnode: " + misclel);
//						System.out.println("Percent for current candidate: " + percent);
						counter = counter + percent;
						wheel[im] = counter;

//						System.out.println("Percentage of hosting for node; " + counter);			
					}
					// Without replacement sorgt fuer einen direkten check ob die knoten ressourcen
					// nicht ausreichen
					if (wheel.length == 1) {

						hostvector[originalPositionOfVNode] = (double) candidates.get(0).getId();
						System.out.println("Treffer! " + hostvector[originalPositionOfVNode]);
						candidates.get(0).setTouched(true);

					}

//					System.out.println(wheel.length);
					for (int is = 1; is < wheel.length; is++) {

						if (rand <= wheel[0]) {

							hostvector[originalPositionOfVNode] = (double) candidates.get(0).getId();
							System.out.println("Treffer! " + hostvector[originalPositionOfVNode]);
							candidates.get(0).setTouched(true);

							break;
						}
						if ((rand > wheel[is - 1]) && rand <= wheel[is]) {

							hostvector[originalPositionOfVNode] = (double) candidates.get(is).getId();
							System.out.println("Treffer! " + hostvector[originalPositionOfVNode]); //
							candidates.get(is).setTouched(true);
							break;
						}

						continue;
					}

					// double misclel = MiscelFunctions.getAr(sn, candidates.get(im));
					// hostprob = misclel/allprob;
					// System.out.println("Allprob " + allprob + " und MiscFunc: " + misclel + "
					// Host: "+ hostprob);
					// if (Math.random() < hostprob) {
					// wenn er hostet:
					// hostvector[originalPositionOfVNode] = (double) candidates.get(im).getId();
					// System.out.println("Treffer! " + hostvector[originalPositionOfVNode]); //
					// candidates.get(im).setTouched(true);
					// break;
					// } else
					// allprob = allprob - misclel;

					// continue;

				}
				break;
			}

			// Schritt4 zweite haelfte
			// System.out.println("Candidate Nodes for: " + vnode.getId() + " " +
			// candidates);
			for (int inn = 0; inn < hostvector.length; inn++) {
//			System.out.println("The host vector an stelle " + inn + ": " + String.valueOf(hostvector[inn]));		
			}
			continue;
		}
		return new Chromosome(c.getLength(), hostvector, statevector);
	}

	private void originalPositionLook(VirtualNetwork vn, Queue<VirtualNode> pq, Map<VirtualNode, Double> pqRanking,
			Map<VirtualNode, Integer> originalPosition, Map<VirtualNode, Integer> mapp) {
		int ist = 0;
		// ALgo arbeitet mit gerichteten Graphen!!
		for (VirtualNode d : pq) {

			Double nr = MiscelFunctions.getAr(vn, d);
			pqRanking.put(d, nr);
			// Well... dont change working code right :)
			originalPosition.put(d, mapp.get(d));
			// Hier ist die Reihenfolge der virtuellen Knoten noch gegeben
			// System.out.println("NodeRank Miscel, ID: " + d.getId() + " and Rank: " + nr);
			// System.out.println("Originale Position von Knoten " + d.getId() + " ist: " +
			// (mapp.get(d)));
			ist++;
			continue;
		}
	}
	// END ALGO1



	private List<SubstrateNode> findFulfillingNodes(AbstractDemand dem) {
		List<SubstrateNode> nodes = new LinkedList<SubstrateNode>();
		for (SubstrateNode n : ns.getSubstrate().getVertices()) {
			if (n.getTouched()) {
				continue;
			} // Falls node touched, i.e. already mapped take next node

			for (AbstractResource res : n)
				if (res.accepts(dem) && res.fulfills(dem)) {
					nodes.add(n);
					break; // Continue with next node.
				}
		}
		return nodes;
	}

	private List<SubstrateNode> findFulfillingNodes2(AbstractDemand dem) {
		List<SubstrateNode> nodes = new LinkedList<SubstrateNode>();
		for (SubstrateNode n : ns.getSubstrate().getVertices()) {

			for (AbstractResource res : n)
				if (res.accepts(dem) && res.fulfills(dem)) {
					nodes.add(n);
					break; // Continue with next node.
				}
		}
		return nodes;
	}

	protected double costOfVN(VirtualNetwork v) {

		double cost = 0;

		for (VirtualNode vn : v.getVertices()) {

			for (AbstractDemand d : vn.get()) {

				if (d instanceof CpuDemand) {

					cost += ((CpuDemand) d).getDemandedCycles();
				}

			}

		}
		for (VirtualLink vl : v.getEdges()) {

			for (AbstractDemand d : vl.get()) {

				if (d instanceof BandwidthDemand) {

					cost += ((BandwidthDemand) d).getDemandedBandwidth();

				}
			}

		}

		return cost;
	}

	@Override // TODO
	protected boolean process(VirtualNetwork v) { // hat als uebergabeparameter das return von getnext()
		double cost = costOfVN(v);
		// Threshhold of 10% for example. If the VN can get embedded with 110% of their
		// minimal cost, then take this solution.
		double threshhold = 1.1 * cost;
		// Wieviele neue generationen
		int iterations = 9; // zu 9
		Chromosome currentBestC = null;
		Double gBest = Double.POSITIVE_INFINITY; // global best
		Double pBest = Double.POSITIVE_INFINITY; // local best
		int popsize = 10; // population size // zu 10
		int counter = v.getVertexCount();
		int var = 0;
		int popCounter = 0;
		this.processedVNetworks++;

		// A Plot for this VN and the fitness of the best-performing member per
		// generation
		// XYSeries series = new XYSeries("Plot of NV: " + v.getLayer());

		VNData datei = new VNData(iterations + 2); // Groesse der CSV datei: 11, davon 10 data und das 11te die
													// minimalen kosten
		datei.writeData(iterations + 1, cost); // minimalen kosten bzw. revenue in 11 der csv datei
		
	//	VNData datei2 = new VNData(iterations + 2);
	//	datei2.writeData(popCounter, cost);
		
		

		Map<VirtualNode, Integer> original = new HashMap<VirtualNode, Integer>();
		// Original hat die originalen Positionen der Knoten und gibt diese an alle
		// anderen methoden weiter
		assignOrder(v, var, original);
		System.out.println("Virtuelle Knoten: " + counter);
		// Initialisierung einer population mit popsize chromosomen und die chr. haben
		// die laenge der VN Knotenanzahl
		Chromosome[] population = new Chromosome[popsize];
		for (int i = 0; i < popsize; i++) {
			population[i] = new Chromosome(counter);
		}
		// The initialization of chromosomes
		for (Chromosome c : population) {

			// c.setZero(); //statevector auf 0 setzen
			c.initializeHosts(); // hostvector auf unendlich setzen
			c = improvedInit(c, this.ns.getSubstrate(), v, original); // Initalisierung der chromosomen bzw. update
			// c.setRandomState(); //Statevector random auf 1 setzen

			// c.printChr();
			// If the first initialization didnt work, the SN ressources are not enough
			for (Double dd : c.getHostvector()) {
				if (dd == Double.POSITIVE_INFINITY) {
					System.out.println(dd);
					System.out.println("NO FEASIBILE SOLUTION! NODE RESSOURCES INSUFFICIENT");
					System.out.println("PROCESSED VNs: ");
					System.out.println(this.processedVNetworks);
					// System.out.println(this.ns.getSubstrate());
					System.out.println("HostingVNs: ");
					System.out.println(this.hostedNetworks);
					return true;
				}
			}
		}

		// The feasibility check and setting the fitness in 1 stage by the 'isFeasibile'
		// method
		for (Chromosome c : population) {
			c.setFeasibile(isFeasibile(this.ns.getSubstrate(), v, c));
			if (c.getFitness() < gBest) {
				System.out.println("besser mit: " + c.getFitness() + " als " + gBest);
				currentBestC = c;
				gBest = c.getFitness();
			}

			c.printChr();
			System.out.println("Feasibile: " + c.getFeasibile());
			System.out.println();
			System.out.println();
		}

		if (currentBestC == null) {

			// series.add(popCounter, 0);
			datei.writeData(0, Double.POSITIVE_INFINITY);
			System.out.println("First Generation was bad.");

		}
		// global best is local best
		// gBest = pBest;
		if (currentBestC != null) {
			System.out.println();
			System.out.println("Gen1 Best Fitness: " + currentBestC.getFitness() + " Equals" + gBest);
			// series.add(popCounter, currentBestC.getFitness());
			currentBestC.printChr();
			// series.add(popCounter, currentBestC.getFitness());
			datei.writeData(0, currentBestC.getFitness());
		}
	//	double sum = 0;
	//	int sumCounter = 0;
	//	for(Chromosome c : population) {
	//		if(c.getFitness() != Double.POSITIVE_INFINITY) {
				
	//			sum += c.getFitness(); 
	//			sumCounter++;
				
				
				
	//		}
			
	//	}
	//	datei2.writeData(popCounter, sum/sumCounter);
		// currentBestC.printChr();
		for (int i = 2; i <= iterations + 1; i++) {

			// Threshhold of 20%: If a good enough solution has been found:break iterations
			if (gBest != Double.POSITIVE_INFINITY && (currentBestC.getFitness()) <= threshhold) {
				broken++;
				break;
			}

			System.out.println();
			System.out.println();
			System.out.println("GENERATION " + i + ":");
			pBest = Double.POSITIVE_INFINITY;
			popCounter++;
			population = generateNewGeneration(population, this.ns.getSubstrate(), v, original);

			// Chromosomen die durch CO oder Mutation entstanden sind wird eine Fitness
			// zugewiesen
			// und das beste Chromosom wird ermittelt
			for (Chromosome chr : population) {

				if (chr.getFitness() == 0.0) {
					System.out.println("Feasibility Check for CO or Mut");
					chr.setFeasibile(isFeasibile(this.ns.getSubstrate(), v, chr));
				}
				// Durch das >= wird der Algo GREEDY jetzt nurnoch >
				if (currentBestC == null || (currentBestC.getFitness() > chr.getFitness() && chr.getFeasibile())) {

					System.out.println("NEW GLOBAL BEST!");
					currentBestC = chr;
					gBest = currentBestC.getFitness();
					System.out.println("NEUE BEST: ");
					currentBestC.printChr();

				}

				if (chr == currentBestC) {
					continue;
				}
//				
				System.out.println("BEFORE new algo1 fitness: " + chr.getFitness());

				// Nochmals alle chromosomen updaten ausser das beste
				// Chromosome clone = doAlgo1(chr, this.ns.getSubstrate(), v, original);
				// clone.setFeasibile(isFeasibile(this.ns.getSubstrate(), v, clone));

				// if (clone.getFitness() < gBest) {

				// currentBestC = clone;
				// gBest = clone.getFitness();

				// }

//				System.out.println("Solo print Gen " + i);
//				chr.printChr();
			}

//			for (Chromosome p : population) {
			// System.out.println("Again just for me: " + p.getFitness());
			//
			// }

			if (currentBestC != null) {
				System.out.println("Gen " + i + " Best mit GBest :" + gBest + " und CBC.fitness: "
						+ currentBestC.getFitness() + " und currentBestC: ");
				currentBestC.printChr();
				// series.add(popCounter, currentBestC.getFitness());
				datei.writeData(popCounter, currentBestC.getFitness());
			}
			
	//		double sum2 = 0;
	//		int sumCounter2 = 0;
	//		for(Chromosome c : population) {
	//			if(c.getFitness() != Double.POSITIVE_INFINITY) {
					
	//				sum2 += c.getFitness(); 
	//				sumCounter2++;
					
					
					
	//			}
				
	//		}
	//		datei2.writeData(popCounter, sum2/sumCounter2);
			
			
		}

		// The local search after the last population is generated
		// for(Chromosome c : population) {
//
//			c = localSearchFirstBest(c, this.ns.getSubstrate(), v, this.ns);
//			if (c.getFitness() < currentBestC.getFitness() && c.getFeasibile()) {
//				
//				currentBestC = c;

//			}

//		}

		if (currentBestC == null || currentBestC.getFitness() == Double.POSITIVE_INFINITY) {

			System.out.println();
//			series.add(popCounter, 0);
			datei.writeData(popCounter, Double.POSITIVE_INFINITY);
			System.out.println("NO FEASIBILE SOLUTION! PROBABLY LINK RESOURCES INSUFFICIENT");

			System.out.println("PROCESSED VNs: ");
			System.out.println(this.processedVNetworks);
			// System.out.println(this.ns.getSubstrate());
			System.out.println("HostingVNs: ");
			System.out.println(this.hostedNetworks);
			return true;

		}

		// XYSeriesCollection data = new XYSeriesCollection(series);
		// JFreeChart chart = ChartFactory.createXYLineChart(

		// "Minimal Fitness Per Generation",
		// "Generation",
		// "minimal fitness value",
		// data,
		// PlotOrientation.VERTICAL,
		// true,
		// true,
		// false

		// );

		// Die naechsten 4 aktivieren und alles davor.
//		ChartPanel chartPanel = new ChartPanel(chart);
//		chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));

//		chartPanel.setDomainZoomable(true);
//		chartPanel.setVisible(true);

		// JPanel pan = new JPanel();

		// pan.add(chartPanel, BorderLayout.CENTER);
		// pan.setLayout(new java.awt.BorderLayout());
		// pan.validate();

		// pan.setVisible(true);

		// Die naechsten 5 aktivieren.
		// JDialog dia = new JDialog();
		// dia.setTitle("Convergence of VN");
		// dia.setSize(450, 300);

		// dia.add(chartPanel);

		// dia.setVisible(true);

		// JFrame frame = new JFrame("Chart");

		// frame.setSize(600, 400);
		// frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		// frame.setVisible(true);
		// frame.getContentPane().add(chartPanel);

		// RefineryUtilities.centerFrameOnScreen();

		// FINAL RETURN

		System.out.println("Definitely the BEST: ");
		// series.add(popCounter, currentBestC.getFitness());
		currentBestC.printChr();

		// public Mapping(AbstractDemand dem, AbstractResource res) {
		// this.demand = dem;
		// demand.register(this);

		// this.resource = res;
		// resource.register(this);
		// }

		// Mapping i.e. Occupying the Node resources of the Chromosome

		Map<VirtualNode, SubstrateNode> finalMap = chrToMap(this.ns.getSubstrate(), v, currentBestC);

		Set<VirtualNode> set = finalMap.keySet();
		Iterator<VirtualNode> virtIt = set.iterator();

		Collection<SubstrateNode> sSet = finalMap.values();
		Iterator<SubstrateNode> subsIt = sSet.iterator();

		// Mapping of the node demands
		while (virtIt.hasNext()) {

			VirtualNode vNode = virtIt.next();
			SubstrateNode sNode = subsIt.next();

			vnm(vNode, sNode);

		}

		// Mapping of the link demands
		AbstractLinkMapping kshort = new kShortestPathLinkMapping(4, false);
		kshort.linkMapping(this.ns.getSubstrate(), v, finalMap);

		// CSV data writer for all minCostValues of all generations
		// final String CSVheader = "Layer, gen1, gen2, gen3, gen4, gen5, gen6, gen7,
		// gen8, gen9, gen10, gen11, gen12, gen13, gen14, gen15, gen16, gen17, gen18,
		// gen19, gen10";
		int counterFile = 0;
		// FileWriter writer = null;

		try {
			// writer = new FileWriter( counterFile + "myData.csv");
			// writer.append(CSVheader);
			// writer.append('\n');

			// writer.append((char)v.getLayer());
			// writer.append(',');

			for (int i = 0; i < datei.getSize(); i++) {
				writer.append(String.valueOf(datei.readData(i)));
				writer.append(',');

			}
			writer.append('\n');

			System.out.println("Writing successful!");
			counterFile++;
		} catch (IOException e) {
			System.out.println("Writing error!");
			e.printStackTrace();
			// } finally {
			// try {
			// writer.flush();
			// writer.close();

			// } catch (IOException e) {
			// System.out.println("Flushing or closing error");
			// e.printStackTrace();
			// }

		}

		// System.out.println("After Process: ");
		this.hostedNetworks++;
		System.out.println("PROCESSED VNs: ");
		System.out.println(this.processedVNetworks);
		// System.out.println(this.ns.getSubstrate());
		System.out.println("HostingVNs: ");
		System.out.println(this.hostedNetworks);
		System.out.println("Cost: " + cost);
		System.out.println("New initializations: " + this.reinit);

		return true;

	}

	private void assignOrder(VirtualNetwork v, int var, Map<VirtualNode, Integer> original) {
		for (VirtualNode vn : v.getVertices()) {

			original.put(vn, var);
			var++;
		}
	}

	// Virtual Node Mapping
	private void vnm(VirtualNode vn, SubstrateNode sn) {

		for (AbstractDemand dem : vn) {

			if (dem instanceof CpuDemand) {
				boolean fulfilled = false;
				for (AbstractResource res : sn) {

					if (res instanceof CpuResource) {
						if (res.accepts(dem) && res.fulfills(dem) && dem.occupy(res)) {
							fulfilled = true;
							break;
						}
						
						if (!fulfilled)
							throw new AssertionError("But we checked before!");
					}
					continue;
				}
			}
			continue;
		}
	}

	// NEW GENERATION

	public Chromosome[] generateNewGeneration(Chromosome[] pop, SubstrateNetwork sn, VirtualNetwork v,
			Map<VirtualNode, Integer> mapp) {
		List<Chromosome> coQ = new LinkedList<Chromosome>(); // Queue for CrossOver and Mutation
		int damnIt = 0; // Counter for unfeasibile chromosomes
		Chromosome best = null; // Lokal bestes chromosom
		double pm = 0.05; // Probability for Mutation //zu 0.05
		double pc = 0.8; // Probability for Crossover
		double fit = Double.POSITIVE_INFINITY;
		Chromosome[] newGeneration = new Chromosome[pop.length];
		Chromosome[] tmpGen = new Chromosome[pop.length];
		double tmpfit = 0;

		// GetBest Chromosome
		for (Chromosome c : pop) {
			System.out.println("This should match the generation before!!" + c.getFitness());
			// isFeasibile(sn, v, c);
			tmpfit = c.getFitness();
			System.out.println("Fitness check: " + tmpfit);
			if (tmpfit < fit && c.getFeasibile()) {
				System.out.println("Besser mit " + c.getFitness());
				fit = tmpfit;
				best = c;
			} else {
				System.out.println("Nicht besser: " + c.getFitness());
				continue;
			}
		}

		// ELITIST
		elitistSelection(pop, v, mapp, best, newGeneration);

		// Alle anderen Chromosomen in eine queue, also ausser dem besten
		for (Chromosome c : pop) {
			//Elitist gone
			if (c == best)
				continue;
			coQ.add(c);
		}
		// System.out.println(tq);

		//

		// Against the comodification error
		Iterator<Chromosome> iterator = coQ.iterator();

		while (iterator.hasNext()) {

			Chromosome c = iterator.next();
			if (c.getFitness() == Double.POSITIVE_INFINITY) {
				System.out.println("Reinitialization");
				this.reinit++;
				// System.out.println("Size before removing: " + coQ.size());
				// c.setZero();
				c = improvedInit(c, sn, v, mapp);
				// c.setRandomState();

//				System.out.println("Is wirklich null hier, bevor es in die naechste Gen kommt" + c.getFitness());
				newGeneration[damnIt] = c;
				// System.out.println("Size after removing: " + coQ.size());
				damnIt++;
				iterator.remove();
			}

		}
		// Siehe Schritt 3 Reinitialization: Alle non feasibile gehen direkt in neachste
		// generation
//		for (Chromosome c : coQ) {
		// Reinitialization Step
		// if (c.getFitness() == Double.POSITIVE_INFINITY) {
		// System.out.println("Size before removing: " + coQ.size());
		// coQ.remove(c);
		// c.setZero();
		// c = doAlgo1(c, sn, v, mapp);
		// c.setRandomState();
		// newGeneration[damnIt] = c;
		// System.out.println("Bad fitness for this chomi! New Initialization and into
		// newGeneration done.");
		// System.out.println("Size before removing: " + coQ.size());
		// coQ.remove(c);
		// System.out.println("Size after removing: " + coQ.size());
		// damnIt++;
		// }
		// }

		System.out.println("How many in coQ: " + coQ.size());

		// Chromosome win1 = tournament(coQ, 2);
		// Chromosome win2 = tournament(coQ, 2);

		// CROSSOVER AND MUTATION
		while (coQ.size() > 1) {

			// Index des zweiten Chromosomes in coQ fuer das Crossover
			int other = generateRandom(1, coQ.size() - 1);
			System.out.println("Das zweite Chromosome fuer eventuellen crossover: " + other);
			System.out.println("Fuer Boundary of: " + coQ.size());

			// Chromosome win1 = tournament(coQ, 2);
			// Chromosome win2 = tournament(coQ, 2);

			// double[] vector1 = win1.getHVector();
			// double[] vector2 = win2.getHVector();

//Die naechsten 2 zeilen wiederherstellen ohne tournament selection
			double[] vector1 = coQ.get(0).getHVector();
			double[] vector2 = coQ.get(other).getHVector();
			// Abbruch da CO mit nur 1 Knoten nicht moeglich
			if (vector1.length == 1) {

				newGeneration[damnIt] = coQ.remove(other);
				damnIt++;
				newGeneration[damnIt] = coQ.remove(0);
				damnIt++;
				continue;
			}
			Queue<Double> GivingNodes = new LinkedList<Double>(); // Nodes die Chr1 hergibt
			Queue<Double> ReceivingNodes = new LinkedList<Double>(); // Nodes die Chr1 bekommt
			ArrayList<Double> CoolNodes1 = new ArrayList<Double>(); // Nodes die Chr1 fuer Infinity einsetzen koennte
																	// ohne
																	// Konflikte
			ArrayList<Double> CoolNodes2 = new ArrayList<Double>(); // Nodes die Chr2 fuer Infinity einsetzen koennte
																	// ohne
																	// Konflikte

			// Original Chromosomes to find possible mutation nodes without crossover
			// Chromosome org1 = coQ.get(0);
			// Chromosome org2 = coQ.get(other);

			// One point crossover happens for 80% of the time

			// if (Math.random() <= pc) {

			// double[] temp = vector1;
			// vector1 = uniform(vector1, vector2);
			// vector2 = uniform(vector2, temp);
			// }

			onePointCO(pc, vector1, vector2, GivingNodes, ReceivingNodes, CoolNodes1, CoolNodes2);

			// List<List<SubstrateNode>> mutationNodes1 = new
			// LinkedList<List<SubstrateNode>>();

			// List<SubstrateNode> possible = new LinkedList<SubstrateNode>();
			// List<SubstrateNode> possible2 = new LinkedList<SubstrateNode>();
			// List<SubstrateNode> gone = new LinkedList<SubstrateNode>();
			// List<SubstrateNode> mutationNodes2 = new LinkedList<SubstrateNode>();
			// original C for finding possible mutation nodes after CO
			Chromosome arg1 = new Chromosome(vector1.length, vector1, coQ.get(0).getSVector());
			Chromosome arg2 = new Chromosome(vector2.length, vector2, coQ.get(other).getSVector());
			System.out.println("THIS NEW");
			arg1.printChr();
			LinkedHashMap<VirtualNode, SubstrateNode> map2 = chrToMap(sn, v, arg2);
			Set<VirtualNode> keys2 = map2.keySet();

			LinkedHashMap<VirtualNode, SubstrateNode> map1 = chrToMap(sn, v, arg1);

			Set<VirtualNode> keys = map1.keySet();

			System.out.println("Key Set: " + keys);
			// Set<VirtualNode> mappedNodes2 = map2.keySet();

			 simpleMutation(pm, vector1, keys);
			 simpleMutation(pm, vector2, keys2);

			// double[] tmp = vector1;
			// System.out.println("ZESS");
			// vector1 = inversionMutReal(pm, vector1, v, sn);
			// vector2 = inversionMutReal(pm, vector2, v, sn);

		//	guidedMutation(pm, vector1, map1, sn);
		//	guidedMutation(pm, vector2, map2, sn);

			// System.out.println("THIS 2 NEW");
			// arg2.printChr();
			// int host2 = 0;

			for (Double d : GivingNodes) {
				System.out.println("GivingNode: " + d);
			}
			for (Double d : ReceivingNodes) {
				System.out.println("ReceivingNode: " + d);
			}

			for (Double d : CoolNodes1) {
				System.out.println("CoolNode1: " + d);
			}
			for (Double d : CoolNodes2) {
				System.out.println("CoolNode2: " + d);
			}

			Chromosome tmp1 = new Chromosome(vector1.length, vector1, coQ.get(0).getSVector());
			Chromosome tmp2 = new Chromosome(vector2.length, vector2, coQ.get(other).getSVector());
			// coQ.get(0).setHostvector(vector1);
			// coQ.get(other).setHostvector(vector2);

			System.out.println("Crossover chromi1: ");
			tmp1.printChr();
			System.out.println("Crossover chromi2: ");
			tmp2.printChr();
			coQ.remove(other);
			coQ.remove(0);
			// feasibile check here

			newGeneration[damnIt] = tmp1;
			damnIt++;
			newGeneration[damnIt] = tmp2;
			damnIt++;
		}

		if (coQ.size() == 1) {

			double[] vector1 = coQ.get(0).getHVector();
			LinkedHashMap<VirtualNode, SubstrateNode> map1 = chrToMap(sn, v, coQ.get(0));

			Set<VirtualNode> keys = map1.keySet();
			simpleMutation(pm, vector1, keys);

			// double[] vector1 = inversionMutReal(pm, last.getHostvector(), v, sn);
			Chromosome arg1 = new Chromosome(vector1.length, vector1, coQ.get(0).getSVector());
			// last.setHostvector(vector1);
			coQ.remove(0);
			newGeneration[damnIt] = arg1;
			// newGeneration[damnIt] = coQ.remove(0);
		}

		 //truncationReplacement(pop, sn, v, tmpGen, newGeneration);

	//	Chromosome[] finalPop = new Chromosome[10];
	//	finalPop = truncReplacement30(pop, newGeneration, sn, v);
		return newGeneration;

	}

	private static double[] inversionMutReal(double pm, double[] vector1, VirtualNetwork vn, SubstrateNetwork sn) {
		double[] result = new double[vector1.length];
//		Iterator<VirtualNode> it = set.iterator();
		VirtualNode testNode = new VirtualNode(2);
		VirtualNode testNode2 = new VirtualNode(2);
		VirtualNode testReal = new VirtualNode(3);
		int vnCounter1 = 0;
		int vnCounter2 = 0;
		VirtualNode myNode1 = null;
		VirtualNode myNode2 = null;
		CpuDemand vnCpu1 = new CpuDemand(0, testNode, "test1");
		CpuDemand vnCpu2 = new CpuDemand(0, testNode2, "test2");
		CpuDemand realDemand = new CpuDemand(0, testReal, "test2");

		for (int i = 0; i < vector1.length; i++) {
			result[i] = vector1[i];
		}

		for (double d : vector1) {

			if (Math.random() <= pm) {

				int trys = 0;
				System.out.print("Try :" + trys + "!");
				while (trys < 5) {

					int pos1 = find(vector1, d);
					int pos2 = generateRandom(0, vector1.length - 1);

					while (pos1 == pos2) {
						pos2 = generateRandom(0, vector1.length - 1);
					}

					System.out.println("swap 1 and 2 :" + pos1 + " " + pos2);
					for (VirtualNode vnode : vn.getVertices()) {
						if (vnCounter1 == find(vector1, d)) {
							myNode1 = vnode;
							for (AbstractDemand dem : myNode1.get()) {

								if (dem instanceof CpuDemand) {
									vnCpu1.setDemandedCycles(((CpuDemand) dem).getDemandedCycles());
									break;
								}
							}
						}
					}
					double d2 = find(vector1, pos2);
					for (VirtualNode vnode : vn.getVertices()) {
						if (vnCounter2 == find(vector1, d2)) {
							myNode2 = vnode;
							for (AbstractDemand dem : myNode2.get()) {

								if (dem instanceof CpuDemand) {
									vnCpu2.setDemandedCycles(((CpuDemand) dem).getDemandedCycles());
									break;
								}
							}
						}
					}

					// Die hostKnoten muessen nur die differenz aktuell hosten koennen
					double realCpu = 0;
					if (vnCpu1.getDemandedCycles() > vnCpu2.getDemandedCycles()) {
						realCpu = vnCpu1.getDemandedCycles() - vnCpu2.getDemandedCycles();
					} else {
						realCpu = vnCpu2.getDemandedCycles() - vnCpu1.getDemandedCycles();
					}
					realDemand.setDemandedCycles(realCpu);

					double id1 = vector1[pos1];
					SubstrateNode candidate1 = null;

					double id2 = vector1[pos2];
					SubstrateNode candidate2 = null;

					for (SubstrateNode sNode : sn.getVertices()) {

						if (sNode.getId() == id2) {

							candidate2 = sNode;
						}
					}
					for (SubstrateNode sNode : sn.getVertices()) {

						if (sNode.getId() == id1) {

							candidate1 = sNode;
						}
					}
					for (AbstractResource res : candidate1)
						if (res.accepts(realDemand) && res.fulfills(realDemand)) {

							for (AbstractResource res2 : candidate2) {

								if (res2.accepts(realDemand) && res2.fulfills(realDemand)) {
									swap(result, pos1, pos2);
									// for(double dd : result) {
									// System.out.println(String.valueOf(dd));
									// }
									mutCheck++;
									return result;

								}
							}

						}
					trys++;
				}
			}
		}
		return result;
	}

	//nicht beachten
	private static int getVnodeAndCPU(AbstractDemand vnCpu1, double[] vector1, VirtualNetwork vn, VirtualNode myNode1,
			double d) {
		int vnCounter1 = 0;
		for (VirtualNode vnode : vn.getVertices()) {
			if (vnCounter1 == find(vector1, d)) {
				myNode1 = vnode;
				for (AbstractDemand dem : myNode1.get()) {

					if (dem instanceof CpuDemand) {
						vnCpu1 = dem;
						break;
					}

				}
				break;
			}
			vnCounter1++;

		}
		return vnCounter1;
	}

	public static int find(double[] a, double target) {

		for (int i = 0; i < a.length; i++) {
			if (a[i] == target)
				return i;
		}
		return -1;
	}

	
	//nicht beachten
	private double[] inversionMut(double[] vector1) {

		double[] result = new double[vector1.length];
		for (int i = 0; i < vector1.length; i++) {
			result[i] = vector1[i];
		}

		int pos1 = generateRandom(0, vector1.length - 1);
		int pos2 = generateRandom(0, vector1.length - 1);

		while (pos1 == pos2) {
			pos2 = generateRandom(0, vector1.length - 1);
		}

		swap(result, pos1, pos2);
		return result;
	}

	private static void swap(double[] vec, int pos1, int pos2) {

		double tmp = vec[pos1];
		vec[pos1] = vec[pos2];
		vec[pos2] = tmp;

	}

	private Chromosome localSearchFirstBest(Chromosome chromi, SubstrateNetwork sn, VirtualNetwork vn,
			NetworkStack ns) {

		long start = System.currentTimeMillis();
		Chromosome fin = new Chromosome(chromi.getLength());
		double actFit = chromi.getFitness();
		// System.out.println("ActFit: " + actFit);
		double[] hosts = chromi.getHostvector();

		// Make CL fuer VNode, exclude all already mapped and calculate fitness. take
		// first which is better.

		int random = generateRandom(0, hosts.length - 1);
		VirtualNode myNode = null;
		int count = 0;

		for (VirtualNode vnode : vn.getVertices()) {
			if (count == random) {
				myNode = vnode;
				break;
			}
			count++;

		}

		List<SubstrateNode> candidates = new LinkedList<SubstrateNode>();

		for (AbstractDemand dem : myNode) {

			if (dem instanceof CpuDemand) {

				candidates = findFulfillingNodesOnly(dem);
				break;
			}
			continue;
		}

		for (SubstrateNode subN : candidates) {

			long id = subN.getId();

			if (check(hosts, id)) {
				continue;
			}
			hosts[count] = id;

			fin.setHostvector(hosts);
			if (isFeasibile(sn, vn, fin) && fin.getFitness() < actFit) {
				timeLocal += System.currentTimeMillis() - start;
				localSuccess++;
				return fin;

			}
			continue;

		}
		// }
		timeLocal += System.currentTimeMillis() - start;
		return chromi;
	}

	private List<SubstrateNode> findFulfillingNodesOnly(AbstractDemand dem) {
		List<SubstrateNode> nodes = new LinkedList<SubstrateNode>();
		int counter = 0;
		for (SubstrateNode n : ns.getSubstrate().getVertices()) {
			if (counter == 10) {
				break;
			}
			if (n.getTouched()) {
				continue;
			} // Falls node touched, i.e. already mapped take next node

			for (AbstractResource res : n)
				if (res.accepts(dem) && res.fulfills(dem)) {
					nodes.add(n);
					counter++;
					break; // Continue with next node.
				}
		}
		return nodes;
	}

	//war mein zuzsaetzlicher crossover operator
	private static double[] anotherCO(double[] vector1, double[] vector2) {

		int cut = generateRandom(1, vector1.length - 1);
		int cut2 = generateRandom(cut, vector1.length - 1) + 1;
		int pos = 0;

		double[] end1 = new double[vector1.length];

		double[] firstSwath = new double[vector1.length];
		double[] secondSwath = new double[vector1.length];

		for (int i = cut; i < cut2; i++) {

			firstSwath[i] = vector1[i];
			end1[i] = vector1[i];

		}
		for (int i = cut; i < cut2; i++) {

			secondSwath[i] = vector2[i];

			// end1[i] = vector1[i];

		}
//		System.out.println("Cut1 :" + cut + " Cut2: " + cut2);
//		System.out.println("Swath1: " + Arrays.toString(firstSwath));
//		System.out.println("Swath2: " + Arrays.toString(secondSwath));

		for (int i = cut; i < cut2; i++) {
			pos = 0;
			if (check(firstSwath, secondSwath[i])) {

				continue;
			}

			double needInsert = secondSwath[i];

			// wenn insert noch nicht in end1 ist, einfuegen.
			if (!check(end1, needInsert)) {

				for (int where = 0; where < end1.length; where++) {

					if (end1[where] == 0) {
						end1[where] = needInsert;
						break;
					}

				}

			}

		}

		for (int in = pos; in < vector1.length; in++) {

			if (end1[in] == 0) {

				for (int is = 0; is < vector2.length; is++) {

					if (!check(end1, vector2[is])) {
						end1[in] = vector2[is];
						//

					}

				}

			}
		}

//		System.out.println("End1: " + Arrays.toString(end1));

		return end1;
	}

	public static boolean check(double[] a, double target) {

		for (int i = 0; i < a.length; i++) {

			if (a[i] == target)
				return true;
		}

		return false;
	}

	private static double[] uniform(double[] vector1, double[] vector2) {

		double[] newUniform = new double[vector1.length];

		for (int i = 0; i < vector1.length; i++) {

			if (Math.random() < 0.1) {

				if (!check(newUniform, vector1[i])) {
					newUniform[i] = vector1[i];
					continue;
				}
			}
			if (!check(newUniform, vector2[i])) {
				newUniform[i] = vector2[i];
				continue;
			}
			if (!check(newUniform, vector1[i])) {
				newUniform[i] = vector1[i];
				continue;
			}
			// Konflikt: Beide wahlmoeglichkeiten sind bereits verbraucht, so nimm den
			// ersten guten aus dem ersten vector.
			for (int counter = 0; counter < vector1.length; counter++) {

				if (!check(newUniform, vector1[counter])) {
					newUniform[i] = vector1[counter];
					break;
				}

			}
		}

		return newUniform;

	}

	private void onePointCO(double pc, double[] vector1, double[] vector2, Queue<Double> GivingNodes,
			Queue<Double> ReceivingNodes, ArrayList<Double> CoolNodes1, ArrayList<Double> CoolNodes2) {
		if (Math.random() <= pc) {

//		System.out.println();
//		System.out.println();
//		System.out.println("CROSSOVER");

			int cut = generateRandom(1, vector1.length - 1);
			System.out.println("Cut-Stelle: " + cut);// Generate cut pointer between 1 and maxlenght-1
			System.out.println();
			for (int i = cut; i < vector1.length; i++) { // Single point crossover
				double tmp = vector1[i];
				GivingNodes.add(tmp);
				vector1[i] = vector2[i];
				ReceivingNodes.add(vector2[i]);
				vector2[i] = tmp;
			}
		}

		// Conflicts back on Infinity
		for (int counter = 0; counter < vector1.length; counter++) {
			for (int counter2 = counter + 1; counter2 < vector1.length; counter2++) {

				// falls ein host 2x vorkommt, wird der zweite wieder auf Infinity gesetzt
				if (vector1[counter] == vector1[counter2] && !(counter == counter2)) {
					vector1[counter2] = Double.POSITIVE_INFINITY;
				}
			}
		}
		for (int counter = 0; counter < vector2.length; counter++) {
			for (int counter2 = counter + 1; counter2 < vector2.length; counter2++) {

				if (vector2[counter] == vector2[counter2] && !(counter == counter2)) {
					vector2[counter2] = Double.POSITIVE_INFINITY;
				}
			}
		}

		// Initialisiere Verfuegbare Knoten als Ersatz fuer Infinity
		for (Double d : GivingNodes) {
			if (ReceivingNodes.contains(d)) {
				continue;
			}
			CoolNodes1.add(d); // Moegliche Knoten um Infinity in Vektor 1 zu ersetzen, CoolNodes2 in Vektor 2
		}
		for (Double d : ReceivingNodes) {
			if (GivingNodes.contains(d)) {
				continue;
			}
			CoolNodes2.add(d);
		}

		for (int i = 0; i < vector1.length; i++) {

			if (vector1[i] == Double.POSITIVE_INFINITY) {
				// TODO
				// Verfuegbare Knoten werden einfach aus der Queue genommen und sind nicht
				// voellig random
				System.out.println("Eintrag in Vektor 1 an Stelle: " + i + " ueberschrieben! Konflikt verhindert.");
				int rand = generateRandom(0, CoolNodes1.size() - 1);
				vector1[i] = CoolNodes1.get(rand);
				CoolNodes1.remove(rand);

			}
		}
		// Infinity entries back to a suitable host-node
		// Idee: Weggetauschte Nodes - Erhaltene Nodes als moegliche Nachfolger fuer
		// Infinity-Eintraege

		for (int i = 0; i < vector2.length; i++) {

			if (vector2[i] == Double.POSITIVE_INFINITY) {
				System.out.println("Eintrag in Vektor 2 an Stelle: " + i + " ueberschrieben! Konflikt verhindert.");
				int rand = generateRandom(0, CoolNodes2.size() - 1);
				vector2[i] = CoolNodes2.get(rand);
				CoolNodes2.remove(rand);
			}
		}
	}

	public Chromosome tournament(List<Chromosome> pop, int k) {

		Chromosome winner = null;
		for (int i = 1; i <= k; i++) {

			Chromosome random = pop.get(generateRandom(0, pop.size()));
			if (winner == null || random.getFitness() < winner.getFitness()) {

				winner = random;
			}
		}
		return winner;
	}

	public Chromosome pmxCO(Chromosome p1, Chromosome p2) {

		Chromosome child = null;
		return child;
	}

	public Chromosome goodMutation(Chromosome chr) {

		Chromosome outcome = null;
		return outcome;

	}

	//nicht beachten
	private void truncationReplacement(Chromosome[] population, SubstrateNetwork sn, VirtualNetwork v,
			Chromosome[] tmpGen, Chromosome[] newGen) {
		Map<Chromosome, Double> sortMe = new HashMap<Chromosome, Double>();
		Map<Chromosome, Double> sorted = new HashMap<Chromosome, Double>();
		for (Chromosome c : newGen) {
			c.setFeasibile(isFeasibile(sn, v, c));
			sortMe.put(c, c.getFitness());

		}
		for (Chromosome c : population) {
			sortMe.put(c, c.getFitness());
		}
		
		//clever, aber leider nimmt er dadurch immer die schlechtesten du honk! sortByValue macht decreasing order, also die schlechtesten
		//chromosome first...
		sorted = MiscelFunctions.sortByValue(sortMe);
		Iterator<Chromosome> chrIt = sorted.keySet().iterator();

		for (int count = 0; count < newGen.length; count++) {
			tmpGen[count] = chrIt.next();

		}
	}
//TODO
	private Chromosome[] truncReplacement30(Chromosome[] oldPopulation, Chromosome[] newPopulation, SubstrateNetwork sn,
			VirtualNetwork v) {

		Chromosome[] result = new Chromosome[oldPopulation.length];
		Map<Chromosome, Double> sortMe = new HashMap<Chromosome, Double>();
		Map<Chromosome, Double> sorted = new HashMap<Chromosome, Double>();
		for (Chromosome c : newPopulation) {
			if (c.getFitness() == 0) {
				c.setFeasibile(isFeasibile(sn, v, c));
			}
			sortMe.put(c, c.getFitness());
		}
		for (Chromosome c : oldPopulation) {
			sortMe.put(c, c.getFitness());
		}
		sorted = MiscelFunctions.sortByValue(sortMe);
		List<Chromosome> list = new LinkedList<Chromosome>();
		List<Chromosome> reverseList = new LinkedList<Chromosome>();
		for(Chromosome c : sorted.keySet()) {
			
			list.add(c);
			
		}
		
		for(int i = list.size()-1; i>= 0; i-- ) {
			
			reverseList.add(list.remove(i));
			
		}
		
		Iterator<Chromosome> chrIt = reverseList.iterator();
		int anzahl = oldPopulation.length/3;
		
		for(int i = 0; i < anzahl; i++) {
			result[i] = chrIt.next();
		}
		int rest = oldPopulation.length-anzahl;
		
		
		for(int c = 0; c < rest; c++) {
			int random = generateRandom(0, newPopulation.length-1);
			
			result[anzahl+c] = newPopulation[random];
			
		}
		return result;
	}

	// Generating a population with tournament selection, mpx crossover, a good
	// mutation mechanism and truncation replacement regarding
	// the new generated chromosomes and list
	//nicht beachten
	public Chromosome[] generateAdvancedGeneration(List<Chromosome> list) {

		Chromosome[] newGen = new Chromosome[list.size()];
		Iterator<Chromosome> iterator = list.iterator();

		while (newGen.length < list.size()) {

			Chromosome parent1 = tournament(list, 2);
			Chromosome parent2 = tournament(list, 2);

			// Chromosome child1 = pmxCO(parent1, parent2);
			// Chromosome child2 = pmxCO(parent2, parent1);

			// Chromosome endChild1 = goodMutation(child1);
			// Chromosome endChild2 = goodMutation(child2);

		}

		return newGen;

	}

	private void elitistSelection(Chromosome[] pop, VirtualNetwork v, Map<VirtualNode, Integer> mapp, Chromosome best,
			Chromosome[] newGeneration) {
		if (best != null) {

			System.out.println("Best found!");
			best.printChr();
			newGeneration[pop.length - 1] = best;
		}
		if (best == null) {

			Chromosome chr = new Chromosome(v.getNodeCount());
			// chr.setZero();
			chr.initializeHosts();
			this.reinit++;
			chr = improvedInit(chr, this.ns.getSubstrate(), v, mapp);
			// chr.setRandomState();
			newGeneration[pop.length - 1] = chr;
		}
	}

	private static List<SubstrateNode> findFulfillingNeighborSet(LinkedHashMap<VirtualNode, SubstrateNode> map1,
			SubstrateNetwork sn, AbstractDemand dem) {
		Set<SubstrateNode> sSet = new HashSet<SubstrateNode>();
		List<SubstrateNode> nope = new LinkedList<SubstrateNode>();

		// in Nodes alle die CPU erfuellen
		List<SubstrateNode> nodes = new LinkedList<SubstrateNode>();
		for (SubstrateNode n : sn.getVertices()) {

			for (AbstractResource res : n)
				if (res.accepts(dem) && res.fulfills(dem)) {
					nodes.add(n);
					break; // Continue with next node.
				}
		}

		// in Sset alle die nachbarn vom bisherigen mapping map1 sind
		for (SubstrateNode n : map1.values()) {

			Collection<SubstrateNode> col = sn.getNeighbors(n);

			for (SubstrateNode snn : col) {

				sSet.add(snn);
			}
		}

		// in nope alle die in nodes, aber nicht in sSet sind
		for (SubstrateNode subN : nodes) {

			if (!(sSet.contains(subN))) {
				nope.add(subN);
			}
		}
		// remove alle knoten in nodes, die in nope sind, also nicht in sSet
		for (SubstrateNode SUBN : nope) {

			nodes.remove(SUBN);

		}

		return nodes;
	}

	// FIXME
	private int guidedMutation(double pm, double[] vector1, LinkedHashMap<VirtualNode, SubstrateNode> map1,
			SubstrateNetwork sn) {
		List<SubstrateNode> possible;
		int host = -1;
		for (VirtualNode virtN : map1.keySet()) {
			host++;
			if (Math.random() < pm) {
				System.out.println("Knoten index: " + String.valueOf(host));
				for (AbstractDemand dem : virtN) {

					if (dem instanceof CpuDemand) {

						possible = findFulfillingNeighborSet(map1, sn, dem); // Hier jetzt 0, 3, 6, 9
						// System.out.println("HEUU");
						// for(SubstrateNode subd : possible) {

						// System.out.println("SubID: " + String.valueOf(subd.getId()));

						// }
						Iterator<SubstrateNode> subsIt = possible.iterator();
						// Iteration ueber die gefundenen Knoten und Elimination aller bereits im
						// hostvector gespeicherten
						// Knoten
						while (subsIt.hasNext()) {

							SubstrateNode snn = subsIt.next();
							// Ueber vektor1 schauen, ob id bereits vergeben
							for (int in = 0; in < vector1.length; in++) {
								if (vector1[in] == snn.getId()) {

									subsIt.remove();
								}
							}
						}

						if (possible.size() > 0) {
							// nimm einen eintrag aus den possible knoten

							// make neighbor list of all mapped SN
							// Set<SubstrateNode> neighbors = neighborSet(map1, sn);
							// System.out.println("SIZE: "+ String.valueOf(possible.size()));
							int rand = generateRandom(0, possible.size() - 1);
							// System.out.println("The random: " + rand);
							double eintrag = possible.get(rand).getId();

							vector1[host] = eintrag;
							possible.remove(rand);

							// System.out.println("Neuer hostvector 1: ");
							// for (Double d : vector1)
							// System.out.println(d);
							mutCheck++;
							break;
						}

						break;
					}
					continue;
				}

			}

			continue;
			// Habe hier also kandidatenknoten fuer virtN

		}
		return host;
	}

	private int simpleMutation2(double pm, double[] vector1, Set<VirtualNode> keys) {
		List<SubstrateNode> possible;
		int host = 0;
		for (VirtualNode virtN : keys) {

			if (Math.random() < pm) {

				for (AbstractDemand dem : virtN) {

					if (dem instanceof CpuDemand) {

						possible = findFulfillingNodes2(dem); // Hier jetzt 0, 3, 6, 9

						Iterator<SubstrateNode> subsIt = possible.iterator();
						// Iteration ueber die gefundenen Knoten und Elimination aller bereits im
						// hostvector gespeicherten
						// Knoten
						while (subsIt.hasNext()) {

							SubstrateNode snn = subsIt.next();
							// Ueber vektor1 schauen, ob id bereits vergeben
							for (int in = 0; in < vector1.length; in++) {
								if (vector1[in] == snn.getId()) {

									subsIt.remove();
								}
							}
						}

						if (possible.size() > 0) {
							// nimm einen eintrag aus den possible knoten

							// make neighbor list of all mapped SN
							// Set<SubstrateNode> neighbors = neighborSet(map1, sn);

							int rand = generateRandom(0, possible.size() - 1);
							System.out.println("The random: " + rand);
							double eintrag = possible.get(rand).getId();

							vector1[host] = eintrag;
							host++;

							possible.remove(rand);

							System.out.println("Neuer hostvector 1: ");
							for (Double d : vector1)
								System.out.println(d);

							break;
						}

						break;
					}
					continue;
				}

			}
			host++;
			continue;
			// Habe hier also kandidatenknoten fuer virtN

		}
		return host;
	}

	private int simpleMutation(double pm, double[] vector1, Set<VirtualNode> keys) {
		List<SubstrateNode> possible;
		int host = 0;
		for (VirtualNode virtN : keys) {

			for (AbstractDemand dem : virtN) {

				if (dem instanceof CpuDemand) {

					possible = findFulfillingNodes2(dem); // Hier jetzt 0, 3, 6, 9

					Iterator<SubstrateNode> subsIt = possible.iterator();

					// Iteration ueber die gefundenen Knoten und Elimination aller bereits im
					// hostvector gespeicherten
					// Knoten
					while (subsIt.hasNext()) {

						SubstrateNode snn = subsIt.next();
						// Ueber vektor1 schauen, ob id bereits vergeben
						for (int in = 0; in < vector1.length; in++) {
							if (vector1[in] == snn.getId()) {

								subsIt.remove();
							}
						}
					}
					System.out.println("The possible vector for 1: " + possible);
					// wenn mutation dann:
					if (Math.random() < pm && possible.size() > 0) {
						// nimm einen eintrag aus den possible knoten
						int rand = generateRandom(0, possible.size() - 1);
						System.out.println("The random: " + rand);
						double eintrag = possible.get(rand).getId();

						vector1[host] = eintrag;
						host++;

						possible.remove(rand);

						System.out.println("Neuer hostvector 1: ");
						for (Double d : vector1)
							System.out.println(d);

						break;
					}
					host++;
					break;
				}
				continue;
			}
			continue;
			// Habe hier also kandidatenknoten fuer virtN

		}
		return host;
	}

	// Random number generation between min and max, both inclusive
	public static int generateRandom(int min, int max) {

		Random random = new Random();
		return random.nextInt(max - min + 1) + min;
	}

	// From Chromosome back to a nodeMapping
	public LinkedHashMap<VirtualNode, SubstrateNode> chrToMap(SubstrateNetwork sn, VirtualNetwork vn, Chromosome c) {
		LinkedHashMap<VirtualNode, SubstrateNode> mapping = new LinkedHashMap<VirtualNode, SubstrateNode>(
				c.getLength());
		Map<VirtualNode, SubstrateNode> mappingArray = new HashMap<>();
		// List<SubstrateNode> subnodes = new LinkedList<SubstrateNode>();
		Collection<VirtualNode> vinodes = new LinkedList<VirtualNode>();
		vinodes = vn.getVertices();
		Iterator<VirtualNode> virtNodes = vinodes.iterator();
		// subnodes = (List<SubstrateNode>) sn.getVertices();
		for (double b : c.getHostvector()) {

			for (SubstrateNode snode : sn.getVertices()) {
				if (snode.getId() == b) {
					VirtualNode next = virtNodes.next();
					mapping.put(next, snode);
					mappingArray.put(next, snode);
					break;
				}
			}
		}
		for (Double d : c.getHostvector()) {

			System.out.println("Chromi: " + d);

		}
		for (VirtualNode d : mapping.keySet()) {

			System.out.println(d.getId() + " " + mapping.get(d).getId());

		}
		return mapping;
	}

//	public SubstrateNetwork copyNetwork(SubstrateNetwork copyFrom) {

	// Neues Netzwerk
//		SubstrateNetwork sn1 = new SubstrateNetwork();
//		for (SubstrateNode sNode : copyFrom.getVertices()) {

	// Fuer jeden knoten einen neuen knoten
//			SubstrateNode SNODE = new SubstrateNode();

	// Fuer jede resource wird sie hinzugefuegt
//			for (AbstractResource res : sNode.get()) {

//				SNODE.add(res.getCopy(sNode, false));

//				for (Mapping m : res.getMappings()) {

//					AbstractDemand d = m.getDemand();
//					
//					for ( AbstractResource RES : SNODE.get()) {
//
//						if ( RES instanceof CpuResource) {
//							
//							d.occupy(RES);	
//						}
//					}
//				}
//			}
//			
//			sn1.addVertex(SNODE);
//		}
//		
//		for (SubstrateLink sLink : copyFrom.getEdges()) {
//				
//			SubstrateLink SLINK = new SubstrateLink();
//			SubstrateNode source = copyFrom.getSource(sLink);
//			SubstrateNode dest = copyFrom.getDest(sLink);
//			
//			SubstrateNode newSource = new SubstrateNode();
//			SubstrateNode newDest = new SubstrateNode();
//			
//			
//			for (SubstrateNode Snode : sn1.getVertices()) {
//				
//				if (Snode.getId() == source.getId()) {
//					
//					newSource = Snode;
//				}
//				if (Snode.getId() == dest.getId()) {
//					
//					newDest = Snode;
//				}	
//			}

//			sn1.addEdge(SLINK, newSource, newDest);
//			
//			for (AbstractResource res : sLink.get()) {
//				
//				
//				
//				SLINK.add(res);
//								
//				for (Mapping m : res.getMappings()) {
//					
//					AbstractDemand d = m.getDemand();
//					
//					for ( AbstractResource RES : SLINK.get()) {
//						
//						if ( RES instanceof BandwidthResource) {
//							
//							d.occupy(RES);	
//						}
//					}
//				}
//			}
//		}	
//		return sn1;
//	}

	// Feasibility-Check: kShortestPath for the Chromosom-Mapping
	public boolean isFeasibile(SubstrateNetwork sn, VirtualNetwork vn, Chromosome c) {

		AbstractLinkMapping kshort = new kShortestPathLinkMappingTest(4, false);
//		AbstractLinkMapping kshort = new kShortestPathLinkMapping(4, true);
		Map<VirtualNode, SubstrateNode> nodeMapping = new LinkedHashMap<VirtualNode, SubstrateNode>();

		Double alreadyMappedLinkCosts = (double) 0;
		BandwidthDemand aDemand;
		// TODO
		SubstrateNetwork testSn = new SubstrateNetwork();

		testSn = sn;

//		System.out.println("TEST SN: " + testSn);

		for (Iterator<SubstrateLink> linkIt = testSn.getEdges().iterator(); linkIt.hasNext();) {
			SubstrateLink currSLink = linkIt.next();
			System.out.println(currSLink);
			AbstractResource res = currSLink.get(BandwidthResource.class);
			// Klar bis hierhin, nur Iteration ueber die SLinks und abgreifen der Bandbreite

			if (res != null) {
				for (Mapping f : res.getMappings()) {
					System.out.println("Ein Mapping! ZUVOR");

					aDemand = (BandwidthDemand) f.getDemand();
					alreadyMappedLinkCosts += aDemand.getDemandedBandwidth();

				}
			}
		}

		nodeMapping = chrToMap(testSn, vn, c);
//		System.out.println(nodeMapping);

		if (kshort.linkMapping(testSn, vn, nodeMapping) && (kshort.processedLinks == kshort.mappedLinks)) {

			Double real = kshort.getCost();
			System.out.println("SHORT - TEST COSTS: " + real);

			// Falls VN nur aus 1 Knoten besteht, ist real = 0
			if (real == 0 && vn.getVertexCount() > 1) {

				System.out.println("REAL is 0!");
				c.setFitness(Double.POSITIVE_INFINITY);

				return false;
			}

			System.out.println("Processed Links: " + kshort.processedLinks);
			System.out.println("Mapped Links " + kshort.mappedLinks + ": Successful!");

			Double linkCost = (double) 0;
			Double allCost = (double) 0;
			Double nodeCost = (double) 0;
			for (VirtualNode vnode : vn.getVertices()) {
				for (AbstractDemand dem : vnode) {
					if (dem instanceof CpuDemand) {
						nodeCost = nodeCost + ((CpuDemand) dem).getDemandedCycles();
					}
					break;
				}

			}
			System.out.println("NodeCost: " + nodeCost);
			BandwidthDemand tmpBwDem;

			// In Berechnung fuer die Kosten: Alle mappings auf den SubstratLinks werden
			// genommen
			for (Iterator<SubstrateLink> tmpSLink = testSn.getEdges().iterator(); tmpSLink.hasNext();) {
				SubstrateLink currSLink = tmpSLink.next();
				System.out.println(currSLink);
				AbstractResource res = currSLink.get(BandwidthResource.class);
				// Klar bis hierhin, nur Iteration ueber die SLinks und abgreifen der Bandbreite

				if (res != null) {
					for (Mapping f : res.getMappings()) {
						System.out.println("Ein Mapping!");

						tmpBwDem = (BandwidthDemand) f.getDemand();
						linkCost += tmpBwDem.getDemandedBandwidth();

					}
				}
			}

			// Falls ein Testmapping erfolgreich war, werden die kosten fuer die bisherigen
			// mappings abgezogen,
			// damit die Kostenrechnung fuer ein einzelnes LinkMapping stimmt
			if (linkCost - alreadyMappedLinkCosts > 0)
				linkCost -= alreadyMappedLinkCosts;

			// LinkCOst wieder ersetzen

			System.out.println("LinkCost: " + real);
			allCost = nodeCost + real;
			c.setFitness(allCost);
			// testSn.clearMappings();
			return true;
		} else
			System.out.println();
		System.out.println("Processed Links: " + kshort.processedLinks);
		System.out.println("Mapped Links " + kshort.mappedLinks + ": Not successful.");
		c.setFitness(Double.POSITIVE_INFINITY);
		// testSn.clearMappings();
		return false;
	}

	@Override
	public List<AbstractAlgorithmStatus> getStati() {

		return null;
	}

	@Override
	protected VirtualNetwork getNext() {
		if (!hasNext())
			return null;

		return curIt.next();
	}

	@SuppressWarnings("unchecked")
	protected boolean hasNext() {
		if (curIt == null || !curIt.hasNext()) {
			if (curNetIt.hasNext()) {
				@SuppressWarnings("unused")
				Network<?, ?, ?> tmp = curNetIt.next();

				// if (tmp instanceof SubstrateNetwork)
				// tmp = curNetIt.next();

				// curIt = this.ns.getVirtuals().iterator();
				curIt = (Iterator<VirtualNetwork>) curNetIt;
				return hasNext();
			} else
				return false;
		} else
			return true;
	}

}
