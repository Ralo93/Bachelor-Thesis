
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
import vnreal.network.NetworkStack;
import vnreal.network.substrate.SubstrateLink;
import vnreal.network.substrate.SubstrateNetwork;
import vnreal.network.substrate.SubstrateNode;
import vnreal.network.virtual.VirtualLink;
import vnreal.network.virtual.VirtualNetwork;
import vnreal.network.virtual.VirtualNode;

public class GeneticAlgoTest {

	public static void main(String[] args) {
		
		
		GeneticAlgoTest test = new GeneticAlgoTest();
		test.algorithmTest();
		
	}

	
	public void algorithmTest() {
		
		
		SubstrateNetwork sn = createSubstrate();
		List<VirtualNetwork> vnw = createVNDemands();
		
		
	//	System.out.println("SN: " + sn);
	//	System.out.println("VN: " + vnw);
		
		NetworkStack stack = new NetworkStack(sn, vnw);
		
	//	System.out.println();

		//System.out.println(stack.getSubstrate());
		//System.out.println(stack.getVirtuals());
		
		IAlgorithm algo = new GeneticAlgo(stack);

		algo.performEvaluation();
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
		cpuRes.setCycles(90.0);
		assertTrue(subsNode1.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode1));

		SubstrateNode subsNode2 = new SubstrateNode();
		idRes = new IdResource(subsNode2, subsNetwork);
		idRes.setId(Integer.toString(2));
		assertTrue(subsNode2.add(idRes));
		cpuRes = new CpuResource(subsNode2);
		cpuRes.setCycles(100.0);
		assertTrue(subsNode2.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode2));

		SubstrateNode subsNode3 = new SubstrateNode();
		idRes = new IdResource(subsNode3, subsNetwork);
		idRes.setId(Integer.toString(3));
		assertTrue(subsNode3.add(idRes));
		cpuRes = new CpuResource(subsNode3);
		cpuRes.setCycles(110.0);
		assertTrue(subsNode3.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode3));

		SubstrateNode subsNode4 = new SubstrateNode();
		idRes = new IdResource(subsNode4, subsNetwork);
		idRes.setId(Integer.toString(4));
		assertTrue(subsNode4.add(idRes));
		cpuRes = new CpuResource(subsNode4);
		cpuRes.setCycles(105.0);
		assertTrue(subsNode4.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode4));

		SubstrateNode subsNode5 = new SubstrateNode();
		idRes = new IdResource(subsNode5, subsNetwork);
		idRes.setId(Integer.toString(5));
		assertTrue(subsNode5.add(idRes));
		cpuRes = new CpuResource(subsNode5);
		cpuRes.setCycles(100.0);
		assertTrue(subsNode5.add(cpuRes));
		assertTrue(subsNetwork.addVertex(subsNode5));

		// Create links
		BandwidthResource bwRes;
		SubstrateLink subsLink12 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink12);
		bwRes.setBandwidth(11.0);
		assertTrue(subsLink12.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink12, subsNode1, subsNode2));

		SubstrateLink subsLink21 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink21);
		bwRes.setBandwidth(11.0);
		assertTrue(subsLink21.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink21, subsNode2, subsNode1));

		SubstrateLink subsLink23 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink23);
		bwRes.setBandwidth(12.0);
		assertTrue(subsLink23.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink23, subsNode2, subsNode3));

		SubstrateLink subsLink32 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink32);
		bwRes.setBandwidth(12.0);
		assertTrue(subsLink32.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink32, subsNode3, subsNode2));

		SubstrateLink subsLink34 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink34);
		bwRes.setBandwidth(9.0);
		assertTrue(subsLink34.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink34, subsNode3, subsNode4));

		SubstrateLink subsLink43 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink43);
		bwRes.setBandwidth(9.0);
		assertTrue(subsLink43.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink43, subsNode4, subsNode3));

		SubstrateLink subsLink45 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink45);
		bwRes.setBandwidth(10.0);
		assertTrue(subsLink45.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink45, subsNode4, subsNode5));

		SubstrateLink subsLink54 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink54);
		bwRes.setBandwidth(10.0);
		assertTrue(subsLink54.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink54, subsNode5, subsNode4));

		SubstrateLink subsLink15 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink15);
		bwRes.setBandwidth(13.0);
		assertTrue(subsLink15.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink15, subsNode1, subsNode5));

		SubstrateLink subsLink51 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink51);
		bwRes.setBandwidth(13.0);
		assertTrue(subsLink51.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink51, subsNode5, subsNode1));

		SubstrateLink subsLink24 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink24);
		bwRes.setBandwidth(15.0);
		assertTrue(subsLink24.add(bwRes));
		assertTrue(subsNetwork.addEdge(subsLink24, subsNode2, subsNode4));

		SubstrateLink subsLink42 = new SubstrateLink();
		bwRes = new BandwidthResource(subsLink42);
		bwRes.setBandwidth(15.0);
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
		vns.add(vn0);

		VirtualNode vA = new VirtualNode(1);
		idDem = new IdDemand(vA);
		idDem.setDemandedId(Integer.toString(1));
		assertTrue(vA.add(idDem));
		cpuDem = new CpuDemand(vA);
		cpuDem.setDemandedCycles(25.0);
		assertTrue(vA.add(cpuDem));
		assertTrue(vn0.addVertex(vA));

		VirtualNode vB = new VirtualNode(1);
		idDem = new IdDemand(vB);
		idDem.setDemandedId(Integer.toString(3));
		assertTrue(vB.add(idDem));
		cpuDem = new CpuDemand(vB);
		cpuDem.setDemandedCycles(27.5);
		assertTrue(vB.add(cpuDem));
		assertTrue(vn0.addVertex(vB));
		
		VirtualNode vY = new VirtualNode(1);
		idDem = new IdDemand(vY);
		idDem.setDemandedId(Integer.toString(6));
		assertTrue(vY.add(idDem));
		cpuDem = new CpuDemand(vY);
		cpuDem.setDemandedCycles(27.5);
		assertTrue(vY.add(cpuDem));
		assertTrue(vn0.addVertex(vY));

		VirtualLink eAB = new VirtualLink(1);
		bwDem = new BandwidthDemand(eAB);
		bwDem.setDemandedBandwidth(9.5);
		assertTrue(eAB.add(bwDem));
		assertTrue(vn0.addEdge(eAB, vA, vB));
		
		VirtualLink eBA = new VirtualLink(1);
		bwDem = new BandwidthDemand(eBA);
		bwDem.setDemandedBandwidth(9.5);
		assertTrue(eBA.add(bwDem));
		assertTrue(vn0.addEdge(eBA, vB, vA));
		
		VirtualLink eYB = new VirtualLink(1);
		bwDem = new BandwidthDemand(eYB);
		bwDem.setDemandedBandwidth(5.5);
		assertTrue(eYB.add(bwDem));
		assertTrue(vn0.addEdge(eYB, vY, vB));

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
		cpuDem.setDemandedCycles(23.5);
		assertTrue(vD.add(cpuDem));
		assertTrue(vn1.addVertex(vD));
		
//		VirtualNode vE = new VirtualNode(2);
//		idDem = new IdDemand(vE);
//		idDem.setDemandedId(Integer.toString(0));
//		assertTrue(vE.add(idDem));
//		cpuDem = new CpuDemand(vE);
//		cpuDem.setDemandedCycles(20.5);
//		assertTrue(vE.add(cpuDem));
//		assertTrue(vn1.addVertex(vE));

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


		return vns;
	}
	
	
	
	
}
