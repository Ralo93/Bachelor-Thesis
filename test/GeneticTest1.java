package vnreal.algorithms;

import static org.junit.Assert.assertTrue;

import java.util.LinkedList;
import java.util.List;

import org.gnu.glpk.GLPK;

import mulavito.algorithms.IAlgorithm;
import vnreal.algorithms.GeneticAlgo;
import vnreal.algorithms.samples.SimpleDijkstraAlgorithm;
import vnreal.constraints.demands.BandwidthDemand;
import vnreal.constraints.demands.CpuDemand;
import vnreal.constraints.demands.IdDemand;
import vnreal.constraints.resources.BandwidthResource;
import vnreal.constraints.resources.CpuResource;
import vnreal.constraints.resources.IdResource;
import vnreal.generators.topologies.WaxmanGenerator;
import vnreal.network.NetworkStack;
import vnreal.network.substrate.SubstrateLink;
import vnreal.network.substrate.SubstrateNetwork;
import vnreal.network.substrate.SubstrateNode;
import vnreal.network.virtual.VirtualLink;
import vnreal.network.virtual.VirtualNetwork;
import vnreal.network.virtual.VirtualNode;

public class GeneticTest1 {


	public static void main(String[] args) {
		
		
	//	WaxmanGenerator generator = new WaxmanGenerator();
		
	//	System.out.println(generator.generateSubstrateNetwork(false));
		
		
		
		
		
		
	GeneticTest1 test = new GeneticTest1();
	test.algorithmTest();
		
	
	
	
	
	
	}

	
	public void algorithmTest() {
		NetworkStack stack = new NetworkStack(createSubstrate(), createVNDemands());

		//System.out.println(stack.getSubstrate());
		//System.out.println(stack.getVirtuals());
		
		IAlgorithm algo = new GeneticAlgo(stack);

		algo.performEvaluation();
		
	//	System.out.println(stack);
		//double f = Double.POSITIVE_INFINITY;
		//System.out.println(f);
	}
	
	
	
	
	
	
	private SubstrateNetwork createSubstrate() {
		SubstrateNetwork subsNetwork = new SubstrateNetwork(true);
		IdResource idRes;
		CpuResource cpuRes;

		SubstrateNode subsNode1 = new SubstrateNode();
		idRes = new IdResource(subsNode1, subsNetwork);
		idRes.setId(Integer.toString(1));
		assertTrue(subsNode1.add(idRes));
		cpuRes = new CpuResource(subsNode1);
		cpuRes.setCycles(100.0);
		assertTrue(subsNode1.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode1));

		SubstrateNode subsNode2 = new SubstrateNode();
		idRes = new IdResource(subsNode2, subsNetwork);
		idRes.setId(Integer.toString(2));
		assertTrue(subsNode2.add(idRes));
		cpuRes = new CpuResource(subsNode2);
		cpuRes.setCycles(60.0);
		assertTrue(subsNode2.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode2));

		SubstrateNode subsNode3 = new SubstrateNode();
		idRes = new IdResource(subsNode3, subsNetwork);
		idRes.setId(Integer.toString(3));
		assertTrue(subsNode3.add(idRes));
		cpuRes = new CpuResource(subsNode3);
		cpuRes.setCycles(30.0);
		assertTrue(subsNode3.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode3));

		SubstrateNode subsNode4 = new SubstrateNode();
		idRes = new IdResource(subsNode4, subsNetwork);
		idRes.setId(Integer.toString(4));
		assertTrue(subsNode4.add(idRes));
		cpuRes = new CpuResource(subsNode4);
		cpuRes.setCycles(20.0);
		assertTrue(subsNode4.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode4));


		// Create links
		BandwidthResource bwRes;
		SubstrateLink subsLink12 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink12);
		bwRes.setBandwidth(20.0);
		assertTrue(subsLink12.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink12, subsNode1, subsNode2));

		SubstrateLink subsLink23 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink23);
		bwRes.setBandwidth(80.0);
		assertTrue(subsLink23.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink23, subsNode2, subsNode3));

		SubstrateLink subsLink31 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink31);
		bwRes.setBandwidth(80.0);
		assertTrue(subsLink31.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink31, subsNode3, subsNode1));

		SubstrateLink subsLink14 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink14);
		bwRes.setBandwidth(70.0);
		assertTrue(subsLink14.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink14, subsNode1, subsNode4));

	
		SubstrateLink subsLink42 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink42);
		bwRes.setBandwidth(70.0);
		assertTrue(subsLink42.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink42, subsNode4, subsNode2));
		

		return subsNetwork;
	}

	private List<VirtualNetwork> createVNDemands() {
		List<VirtualNetwork> vns = new LinkedList<VirtualNetwork>();
		IdDemand idDem;
		CpuDemand cpuDem;
		BandwidthDemand bwDem;

		// Virtual network 0
		VirtualNetwork vn0 = new VirtualNetwork(1);
	//	vns.add(vn0);

		VirtualNode vJ = new VirtualNode(1);
		idDem = new IdDemand(vJ);
		idDem.setDemandedId(Integer.toString(1));
		assertTrue(vJ.add(idDem));
		cpuDem = new CpuDemand(vJ);
		cpuDem.setDemandedCycles(15.0);
		assertTrue(vJ.add(cpuDem));
		assertTrue(vn0.addVertex(vJ));

		VirtualNode vB = new VirtualNode(1);
		idDem = new IdDemand(vB);
		idDem.setDemandedId(Integer.toString(3));
		assertTrue(vB.add(idDem));
		cpuDem = new CpuDemand(vB);
		cpuDem.setDemandedCycles(10.0);
		assertTrue(vB.add(cpuDem));
		assertTrue(vn0.addVertex(vB));
		
		
		VirtualLink eAB = new VirtualLink(1);
		bwDem = new BandwidthDemand(eAB);
		bwDem.setDemandedBandwidth(10.0);
		assertTrue(eAB.add(bwDem));
		assertTrue(vn0.addEdge(eAB, vJ, vB));
		
		VirtualLink eBA = new VirtualLink(1);
		bwDem = new BandwidthDemand(eBA);
		bwDem.setDemandedBandwidth(5.0);
		assertTrue(eBA.add(bwDem));
		assertTrue(vn0.addEdge(eBA, vB, vJ));
		
		// Virtual network 1
				VirtualNetwork vn1 = new VirtualNetwork(2);
				vns.add(vn1);

				VirtualNode vC = new VirtualNode(2);
				idDem = new IdDemand(vC);
				idDem.setDemandedId(Integer.toString(5));
				assertTrue(vC.add(idDem));
				cpuDem = new CpuDemand(vC);
				cpuDem.setDemandedCycles(25.0);
				assertTrue(vC.add(cpuDem));
				assertTrue(vn1.addVertex(vC));

				VirtualNode vD = new VirtualNode(2);
				idDem = new IdDemand(vD);
				idDem.setDemandedId(Integer.toString(2));
				assertTrue(vD.add(idDem));
				cpuDem = new CpuDemand(vD);
				cpuDem.setDemandedCycles(20.0);
				assertTrue(vD.add(cpuDem));
				assertTrue(vn1.addVertex(vD));
				
				VirtualNode vE = new VirtualNode(2);
				idDem = new IdDemand(vE);
				idDem.setDemandedId(Integer.toString(0));
				assertTrue(vE.add(idDem));
				cpuDem = new CpuDemand(vE);
				cpuDem.setDemandedCycles(21.0);
				assertTrue(vE.add(cpuDem));
				assertTrue(vn1.addVertex(vE));

				VirtualLink eCD = new VirtualLink(2);
				bwDem = new BandwidthDemand(eCD);
				bwDem.setDemandedBandwidth(9.0);
				assertTrue(eCD.add(bwDem));
				assertTrue(vn1.addEdge(eCD, vC, vD));
				
				VirtualLink eDC = new VirtualLink(2);
				bwDem = new BandwidthDemand(eDC);
				bwDem.setDemandedBandwidth(9.0);
				assertTrue(eDC.add(bwDem));
				assertTrue(vn1.addEdge(eDC, vD, vC));
				
				//VN 2
				VirtualNetwork vn2 = new VirtualNetwork(3);
			//	vns.add(vn2);

				VirtualNode vP = new VirtualNode(3);
				idDem = new IdDemand(vP);
				idDem.setDemandedId(Integer.toString(1));
				assertTrue(vP.add(idDem));
				cpuDem = new CpuDemand(vP);
				cpuDem.setDemandedCycles(3.5);
				assertTrue(vP.add(cpuDem));
				assertTrue(vn2.addVertex(vP));

				VirtualNode vK = new VirtualNode(3);
				idDem = new IdDemand(vK);
				idDem.setDemandedId(Integer.toString(3));
				assertTrue(vK.add(idDem));
				cpuDem = new CpuDemand(vK);
				cpuDem.setDemandedCycles(2.0);
				assertTrue(vK.add(cpuDem));
				assertTrue(vn2.addVertex(vK));
				
				
				VirtualLink ePK = new VirtualLink(3);
				bwDem = new BandwidthDemand(ePK);
				bwDem.setDemandedBandwidth(12.0);
				assertTrue(ePK.add(bwDem));
				assertTrue(vn2.addEdge(ePK, vP, vK));
		

		

		return vns;
	}
	
	
	
	
}

