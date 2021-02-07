/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2010 Network Security Lab, University of Washington, Seattle.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "ns3/config-store-module.h"
#include "ns3/wifi-module.h"
#include "ns3/energy-module.h"
#include "ns3/internet-module.h"
#include "ns3/reliability-module.h"
#include "ns3/applications-module.h"
#include "ns3/internet-module.h"
#include "ns3/olsr-helper.h"
#include "ns3/flow-monitor-helper.h"
#include "ns3/flow-monitor-module.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("ReliabilityExample");

std::vector < std::string > split(std::string const & str, const char delim);
double compute_soh(std::vector<double> Tc);
double compute_mttf(std::vector<double> Tc);

/**
 * \param socket Pointer to socket.
 *
 * Packet receiving sink.
 */
void
ReceivePacket (Ptr<Socket> socket)
{
  Ptr<Packet> packet;
  Address from;
  while ((packet = socket->RecvFrom (from)))
    {
      if (packet->GetSize () > 0)
        {
          InetSocketAddress iaddr = InetSocketAddress::ConvertFrom (from);

          std::cout << " Node " << socket->GetNode()->GetId()
                    << " received one packet from Socket: " << iaddr.GetIpv4 ()
                    << " port: " << iaddr.GetPort ()
                    << " size: " << packet->GetSize ()
                    << " at time = " << Simulator::Now ().GetSeconds ()
                    << std::endl;
        }
    }
}

static void SendPacket (Ptr<Socket> socket, uint32_t pktSize, Time pktInterval )
{
  socket->Send (Create<Packet> (pktSize));
  Simulator::Schedule (pktInterval, &SendPacket, socket, pktSize, pktInterval);
}

/// Trace function for remaining energy at node.
void
RemainingEnergy (double oldValue, double remainingEnergy)
{
  NS_LOG_UNCOND (Simulator::Now ().GetSeconds ()
                 << "s Current remaining energy = " << remainingEnergy << "J");
}

/// Trace function for average power consumption at node.
void
AveragePower (double oldValue, double totalEnergy)
{
  NS_LOG_UNCOND (Simulator::Now ().GetSeconds ()
                 << "s Average Power = " << totalEnergy/Simulator::Now ().GetSeconds () << "W");
}

/// Trace function for total energy consumption at node.
void
TotalEnergy (double oldValue, double totalEnergy)
{
  NS_LOG_UNCOND (Simulator::Now ().GetSeconds ()
                 << "s Total energy consumption = " << totalEnergy << "J");
}

/// Trace function for the power harvested by the energy harvester.
void
HarvestedPower (double oldValue, double harvestedPower)
{
  std::cout << Simulator::Now ().GetSeconds ()
            << "s Current harvested power = " << harvestedPower << " W" << std::endl;
}

/// Trace function for the total energy harvested by the node.
void
TotalEnergyHarvested (double oldValue, double TotalEnergyHarvested)
{
  std::cout << Simulator::Now ().GetSeconds ()
            << "s Total energy harvested by harvester = "
            << TotalEnergyHarvested << " J" << std::endl;
}

void
PrintAveragePower (std::string pwrFile, EnergySourceContainer sources, int simSpeed)
{
  std::fstream results;
  results.open (pwrFile, std::ios::app ); 
	
  uint32_t nSources = sources.GetN ();
  Ptr<BasicEnergySource> basicSourcePtr;
  Ptr<DeviceEnergyModel> basicRadioModelPtr;
  for (uint8_t i = 0; i < nSources; ++i)
  {
    basicSourcePtr = DynamicCast<BasicEnergySource> (sources.Get (i));
    basicRadioModelPtr = basicSourcePtr->FindDeviceEnergyModels ("ns3::CpuEnergyModel").Get (0);
    // Ptr<EnergySource> sourcei = sources.Get(i);
    std::cout << " At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<basicSourcePtr->GetNode()->GetId();
    std::cout << " Average Power = " << basicRadioModelPtr->GetTotalEnergyConsumption()/(double(simSpeed)*Simulator::Now ().GetSeconds ())<<std::endl;
    results << basicRadioModelPtr->GetTotalEnergyConsumption()/(double(simSpeed)*Simulator::Now ().GetSeconds ()) << "\n";
  }
  results.close ();
}

void
PrintEnergy (std::string engFile, EnergySourceContainer sources, double secondsPeriod)
{
  std::fstream results;
  results.open (engFile, std::ios::app ); 
	
  uint32_t nSources = sources.GetN ();
  for (uint8_t i = 0; i < nSources; ++i)
  {
    Ptr<EnergySource> sourcei = sources.Get(i);
    // std::cout << " At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<sourcei->GetNode()->GetId();
    // std::cout << " Remaining Capacity = " << sourcei->GetRemainingEnergy()/(3.6*3.3)<<std::endl;
    results << sourcei->GetRemainingEnergy()/(3.6*3.3) << " ";
  }
  results << "\n";
  results.close ();
  if (!Simulator::IsFinished ())
  {
    Simulator::Schedule (Seconds (secondsPeriod), &PrintEnergy,engFile,sources,secondsPeriod);
  }
}

void
PrintReliability (std::string relFile, NodeContainer nodes, double secondsPeriod)
{
  std::fstream results;
  results.open (relFile, std::ios::app ); 

  uint32_t nDevices = nodes.GetN ();
  for (uint8_t i = 0; i < nDevices; ++i)
  {
    Ptr<Node> nodei = nodes.Get(i);
    // std::cout << " At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<nodei->GetId();
    // std::cout << " Reliability = " << nodei->GetObject<ReliabilityModel>()->GetReliability()<<std::endl;
    results << nodei->GetObject<ReliabilityModel>()->GetReliability() << " ";
  }
  results << "\n";
  results.close ();
  if (!Simulator::IsFinished ())
  {
    Simulator::Schedule (Seconds (secondsPeriod), &PrintReliability,relFile,nodes,secondsPeriod);
  }
}

void
PrintTime ()
{
    std::cout << " Simulation Time: "<< Simulator::Now().GetSeconds()<<std::endl;
  if (!Simulator::IsFinished ())
  {
    Simulator::Schedule (Seconds (100), &PrintTime);
  }
}
// void
// PrintEnergy (std::ofstream& results, EnergySourceContainer sources, double secondsPeriod)
// {

//   uint32_t nSources = sources.GetN ();
//   for (uint8_t i = 0; i < nSources; ++i)
//   {
//     Ptr<EnergySource> sourcei = sources.Get(i);
//     std::cout << " At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<sourcei->GetNode()->GetId();
//     std::cout << " Remaining Energy = " << sourcei->GetRemainingEnergy()<<std::endl;
//     results << sourcei->GetRemainingEnergy() << " ";
//   }
//   results << "\n";
//   // if (!Simulator::IsFinished ())
//   // {
//   //   Simulator::Schedule (Seconds (secondsPeriod), &PrintEnergy,results,sources,secondsPeriod);
//   // }
// }

// void
// PrintReliability (std::ofstream& results, NodeContainer nodes, double secondsPeriod)
// {
//   uint32_t nDevices = nodes.GetN ();
//   for (uint8_t i = 0; i < nDevices; ++i)
//   {
//     Ptr<Node> nodei = nodes.Get(i);
//     std::cout << " At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<nodei->GetId();
//     std::cout << " Reliability = " << nodei->GetObject<ReliabilityModel>()->GetReliability()<<std::endl;
//     results << nodei->GetObject<ReliabilityModel>()->GetReliability() << " ";
//   }
//   results << "\n";
//   // if (!Simulator::IsFinished ())
//   // {
//   //   Simulator::Schedule (Seconds (secondsPeriod), &PrintReliability,results,nodes,secondsPeriod);
//   // }
// }


void
PrintInfo (Ptr<Node> node)
{

  std::cout << " At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<node->GetId();
  std::cout << " CPU Power = " << node->GetObject<PowerModel>()->GetPower();
  std::cout << " Performance = " << node->GetObject<PerformanceModel>()->GetThroughput();
  std::cout << " Temperature = " << node->GetObject<TemperatureModel>()->GetTemperature()<<std::endl;
  std::cout << " Reliability = " << node->GetObject<ReliabilityModel>()->GetReliability()<<std::endl;

  if (!Simulator::IsFinished ())
  {
    Simulator::Schedule (Seconds (0.5), &PrintInfo,node);
  }
}

// Inputs
std::string srFile = "sr_15_OPT_wo.txt";     // Sensor location file
std::string gwFile = "gw_15_OPT_wo.txt";      // Gateway location file
std::string flFile = "fl_15_OPT_wo.txt";      // File of transmission distance at each node
std::string ptxFile = "ptx_15_OPT_wo.txt";    // File of transmission power at each node
std::string tempFile = "temp_15_OPT_wo.txt";  // Temperature traces at each sensor location
std::string rcFile = "rc_15_OPT_wo.txt";     // Recharge power at each sensor node

// Outputs
std::string relFile = "rel";     // Reliability for all nodes output
std::string engFile = "eng";     // Remaining energy for all nodes output

std::string sohFile = "soh";     // SoH for all nodes output
std::string mttfFile = "mttf";     // MTTF for all nodes output

std::string pwrFile = "pwr";     // Average power for all nodes output


int
main (int argc, char *argv[])
{
  LogComponentEnable ("BasicEnergyHarvester", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("CpuEnergyModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("StatePowerModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("PowerModel", LOG_LEVEL_DEBUG);
  //LogComponentEnable ("TemperatureSimpleModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("ReliabilityTDDBModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("AppPowerModel", LOG_LEVEL_DEBUG);
  
  uint32_t nDevices = 0;
  uint32_t nGateways = 0;
  std::string phyMode ("DsssRate1Mbps");
  double Prss = -80;            // dBm
  uint32_t packetSize = 1000;   // bytes
  uint32_t originalPacketSize = 100;   // bytes
  std::string dataRate = "200Bps";                  /* Application layer datarate. */

  double bw = 2000.0;             // B/s
  // uint32_t dataSize = 10000;    // bytes for reliability helper
  //std::string dataRate = "2000Bps";                  

  // Energy parameters
  double Es = 0.04;             // J, energy for one sensing sample
  double P0 = 0.01;             // W, ambient power dissipation
  // double Pto = 0.22;            // W, ambient wifi transmission power
  // double alpha = 3.5;           // Path exponent
  // double beta = 0.0000001;      // Path linear coefficient
  double Prx = 0.1;             // W, wifi receiption power
  double BattJ = 1188000;         // J, battery capacity - 10,000mAh - 118,800 J

  // Simulation parameters
  double packetInterval = 5;    // seconds
  double startTime = 0.0;       // seconds
  uint32_t simSpeed = 3600;  // i.e. 600 for 1sec simulation == 10min real time 
  double simulationDays = 12*30.0;  // days
  double simulationSeconds = simulationDays*24*60*60/simSpeed;
  // uint32_t port = 5000;         // port number that receiver is listening on
  bool verbose = false;
  bool tracing = true;
  bool flow_monitor = true;
  // Energy Harvester variables
  double harvestingUpdateInterval = 0.05;  // seconds
  /*
   * This is a magic number used to set the transmit power, based on other
   * configuration.
   */
  double offset = 81;
  std::string txtName = "15_OPT_wo";  

  CommandLine cmd;
  cmd.AddValue ("phyMode", "Wifi Phy mode", phyMode);
  cmd.AddValue ("Prss", "Intended primary RSS (dBm)", Prss);
  cmd.AddValue ("packetSize", "size of application packet sent", packetSize);
  cmd.AddValue ("startTime", "Simulation start time", startTime);
  cmd.AddValue ("txtName", "Name of input text files", txtName);
  cmd.Parse (argc, argv);

  srFile = "sr_" + txtName + ".txt";     // Sensor location file
  gwFile = "gw_" + txtName + ".txt";      // Gateway location file
  flFile = "fl_" + txtName + ".txt";      // File of transmission distance at each node
  ptxFile = "ptx_" + txtName + ".txt";    // File of transmission power at each node
  tempFile = "temp_" + txtName + ".txt"; // Temperature traces at each sensor location
  rcFile = "rc_" + txtName + ".txt";     // Recharge power at each sensor node
  
  relFile = "rel_" + txtName + ".txt";     // Reliability for all nodes output
  engFile = "eng_" + txtName + ".txt";     // Remaining energy for all nodes output
  sohFile = "soh_" + txtName + ".txt";     // SoH for all nodes output
  mttfFile = "mttf_" + txtName + ".txt";     // MTTF  for all nodes output
  pwrFile = "pwr_" + txtName + ".txt";     // Average power for all nodes output

  std::ofstream reliabilityResults;
  reliabilityResults.open(relFile);
  std::ofstream energyResults;
  energyResults.open(engFile);
  std::ofstream powerResults;
  powerResults.open(pwrFile);


  reliabilityResults.close ();
  energyResults.close ();
  powerResults.close ();

  // disable fragmentation for frames below 2200 bytes
  Config::SetDefault ("ns3::WifiRemoteStationManager::FragmentationThreshold",
                      StringValue ("2200"));
  // turn off RTS/CTS for frames below 2200 bytes
  Config::SetDefault ("ns3::WifiRemoteStationManager::RtsCtsThreshold",
                      StringValue ("2200"));
  // Fix non-unicast data rate to be the same as that of unicast
  Config::SetDefault ("ns3::WifiRemoteStationManager::NonUnicastMode",
                      StringValue (phyMode));


  ///////////////////////////
  // Create sensor devices //
  ///////////////////////////

  NodeContainer sensorDevices;
  MobilityHelper mobilityEd;
  Ptr<ListPositionAllocator> positionAllocEd = CreateObject<ListPositionAllocator> ();

  // Read end nodes' locations from text file
  std::ifstream EdLocationFile(srFile);
  std::vector<int> SensorFlag;       // Whether the node is a sensor node (1) or a relay node (0)
  if (EdLocationFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing sensor location file.");
    std::string line;
    while (std::getline(EdLocationFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > coordinates = split(line, ' ');
            double x = atof(coordinates.at(0).c_str());
            double y = atof(coordinates.at(1).c_str());
            int Sensor = atof(coordinates.at(2).c_str());
            positionAllocEd->Add (Vector (x, y, 0.0) );
            SensorFlag.push_back(Sensor);
            nDevices ++;
        }
    }
    mobilityEd.SetPositionAllocator (positionAllocEd);
    mobilityEd.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << srFile);
    return -1;
  }

  sensorDevices.Create (nDevices);
  mobilityEd.Install (sensorDevices);

  ////////////////
  // Create GWs //
  ////////////////

  NodeContainer gateways;

  MobilityHelper mobilityGw;
  Ptr<ListPositionAllocator> positionAllocGw = CreateObject<ListPositionAllocator> ();

  // Read gateway locations from text file
  std::ifstream GwLocationFile(gwFile);
  if (GwLocationFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing gw device location file.");
    std::string line;
    while (std::getline(GwLocationFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > coordinates = split(line, ' ');
            double x = atof(coordinates.at(0).c_str());
            double y = atof(coordinates.at(1).c_str());
            positionAllocGw->Add (Vector (x, y, 15.0) );
            nGateways ++;
        }
    }
    mobilityGw.SetPositionAllocator (positionAllocGw);
    mobilityGw.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << gwFile);
    return -1;
  }

  gateways.Create (nGateways);
  mobilityGw.Install (gateways);

  NodeContainer allNodes = NodeContainer(sensorDevices, gateways);

  //////////////////////
  // Read flow matrix //
  //////////////////////
  std::ifstream FLFile(flFile);
  std::vector< std::vector<double> > fij(nDevices); // flow matrix
  if (FLFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing flow file.");
    std::string line;
    int srcIdx = 0;
    while (std::getline(FLFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > fijLine = split(line, ' ');
            
            for (std::vector < std::string >::iterator it = fijLine.begin();
              it != fijLine.end(); ++it)
            {
              double fijVal = atof(it->c_str());
              fij[srcIdx].push_back(fijVal);
            }
            srcIdx ++;
        }
    }
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << flFile);
    return -1;
  }


  /////////////////////////////////////
  // Read ambient transmission power //
  /////////////////////////////////////
  std::ifstream PtxFile(ptxFile);
  std::vector<double> ptxo; // flow matrix
  if (PtxFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing distance file.");
    std::string line;
    while (std::getline(PtxFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > data = split(line, ' ');
            ptxo.push_back(atof(data.at(0).c_str()));
        }
    }
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << ptxFile);
    return -1;
  }


 //////////////////////
  // Read temperature matrix //
  //////////////////////
  std::ifstream TempFile(tempFile);
  std::vector< std::vector<double> > tempMatrix(nDevices); // temperature matrix
  if (TempFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing temperature file.");
    std::string line;
    int srcIdx = 0;
    while (std::getline(TempFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > tempLine = split(line, ' ');
            
            for (std::vector < std::string >::iterator it = tempLine.begin();
              it != tempLine.end(); ++it)
            {
              double tempVal = atof(it->c_str());
              tempMatrix[srcIdx].push_back(tempVal);
            }
            srcIdx ++;
        }
    }
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << tempFile);
    return -1;
  }


 //////////////////////
  // Read recharge power matrix //
  //////////////////////
  std::ifstream RcFile(rcFile);
  std::vector< std::vector<double> > rcMatrix(nDevices); // recharge matrix
  if (RcFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing recharge file.");
    std::string line;
    int srcIdx = 0;
    while (std::getline(RcFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > rcLine = split(line, ' ');
            
            for (std::vector < std::string >::iterator it = rcLine.begin();
              it != rcLine.end(); ++it)
            {
              double rcVal = atof(it->c_str())/(1.23);
              rcMatrix[srcIdx].push_back(rcVal);
            }
            srcIdx ++;
        }
    }
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << rcFile);
    return -1;
  }


  ////////////////////////
  // Wifi Configuration //
  ////////////////////////
  WifiHelper wifi;
  if (verbose)
  {
    wifi.EnableLogComponents ();  // Turn on all Wifi logging
  }

  YansWifiPhyHelper wifiPhy =  YansWifiPhyHelper::Default ();
  // set it to zero; otherwise, gain will be added
  wifiPhy.Set ("RxGain", DoubleValue (-10) );
  wifiPhy.Set ("TxGain", DoubleValue (offset + Prss));
  wifiPhy.Set ("CcaMode1Threshold", DoubleValue (0.0));
  // ns-3 supports RadioTap and Prism tracing extensions for 802.11b
  wifiPhy.SetPcapDataLinkType (YansWifiPhyHelper::DLT_IEEE802_11_RADIO);

  YansWifiChannelHelper wifiChannel;
  wifiChannel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
  wifiChannel.AddPropagationLoss ("ns3::FriisPropagationLossModel");
  wifiPhy.SetChannel (wifiChannel.Create ());

  // Disable rate control
  wifi.SetStandard (WIFI_PHY_STANDARD_80211b);
  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
                                "DataMode",StringValue (phyMode),
                                "ControlMode",StringValue (phyMode));
  // Set it to adhoc mode
  WifiMacHelper wifiMac;
  wifiMac.SetType ("ns3::AdhocWifiMac");

  /** install PHY + MAC **/
  NetDeviceContainer devicesNet = wifi.Install (wifiPhy, wifiMac, sensorDevices);
  NetDeviceContainer gatewaysNet = wifi.Install (wifiPhy, wifiMac, gateways);
  NetDeviceContainer allNodesNet = NetDeviceContainer(devicesNet, gatewaysNet);

  /** Internet stack **/
  InternetStackHelper internet;
  internet.Install (allNodes);

  Ipv4AddressHelper ipv4;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer allNodesInterface = ipv4.Assign (allNodesNet);

  // Configure the static routing from the flow matrix
  Ipv4StaticRoutingHelper staticRouting;
  // Ipv4ListRoutingHelper list;
  // list.Add (staticRouting, 0);

  internet.SetRoutingHelper (staticRouting); // has effect on the next Install ()

  /** Socket-Source/Sink **/
  /***************************************************************************/
  ApplicationContainer sourceApplications, sinkApplications;
  uint32_t portNumber = 9;

  for (uint8_t i = 0; i < nDevices; ++i)
  {
    for (uint8_t j = 0; j < nDevices + nGateways; ++j)
    {
      if (fij[i][j] < 0.1) // If no flow is assigned, skip this connection
        continue;

      Ptr<Node>srNode = sensorDevices.Get(i);
      Ptr<Ipv4> ipv4Node = srNode->GetObject<Ipv4> ();
      Ptr<Ipv4StaticRouting> staticRoutingNode = staticRouting.GetStaticRouting (ipv4Node);
      Ipv4Address sinkAddr = allNodesInterface.GetAddress(j);
      Ipv4Address nextHopAddr = allNodesInterface.GetAddress(j);
      staticRoutingNode->AddHostRouteTo (sinkAddr, nextHopAddr, 1);
      // sinkAddr.Print(std::cout);
      // nextHopAddr.Print(std::cout);
      auto ipv4 = allNodes.Get(j)->GetObject<Ipv4> (); //sink
      const auto address = ipv4->GetAddress (1, 0).GetLocal ();
      InetSocketAddress sinkSocket (address, portNumber++);
      OnOffHelper onOffHelper ("ns3::UdpSocketFactory", sinkSocket);
      onOffHelper.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
      onOffHelper.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
      onOffHelper.SetAttribute ("DataRate", DataRateValue (8*(6.25*fij[i][j])*packetSize/(1*20*5)));
      onOffHelper.SetAttribute ("PacketSize", UintegerValue (packetSize)); //bytes
      sourceApplications.Add (onOffHelper.Install (allNodes.Get (i)));
      PacketSinkHelper packetSinkHelper ("ns3::UdpSocketFactory", sinkSocket);
      sinkApplications.Add (packetSinkHelper.Install (allNodes.Get (j)));
      
    }
  }

  // sink = StaticCast<PacketSink> (sinkApplications.Get (0));
  sinkApplications.Start (Seconds (0.0));
  sourceApplications.Start (Seconds (0.0));

  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();


  // TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  // // Configure receiving node - the only sink
  // Ipv4Address sinkAddr = allNodesInterface.GetAddress(nDevices+nGateways-1);
  // Ptr<Socket> recvSocket = Socket::CreateSocket (allNodes.Get(nDevices+nGateways-1), tid);
  // InetSocketAddress local = InetSocketAddress (sinkAddr, port);
  // recvSocket->Bind (local);
  // //recvSocket->SetRecvCallback (MakeCallback (&ReceivePacket));

  // // Configure source node
  // for (uint8_t i = 0; i < nDevices; ++i)
  // {
  //   if (SensorFlag[i] > 0)
  //   {
  //     Ptr<Socket> source = Socket::CreateSocket (sensorDevices.Get(i), tid);
  //     InetSocketAddress remote = InetSocketAddress (sinkAddr, port);
  //     source->Connect (remote);

  //     // Schedule SendPacket
  //     Simulator::ScheduleWithContext(source->GetNode()->GetId(),
  //                                    Seconds (startTime), &SendPacket,
  //                                    source, packetSize,
  //                                    Seconds(packetInterval));
  //   }
  // }

  ////////////////////////////////////////
  // Power, Temperature and Reliability //
  ////////////////////////////////////////
  /** Energy Model **/
  /***************************************************************************/
  /* Energy source */
  BasicEnergySourceHelper basicSourceHelper;
  // Configure energy source
  basicSourceHelper.Set ("BasicEnergySourceInitialEnergyJ", DoubleValue (BattJ));
  basicSourceHelper.Set ("BasicEnergySupplyVoltageV", DoubleValue (3.3));
  // Install source
  EnergySourceContainer energySources = basicSourceHelper.Install (sensorDevices);




  /* Reliability stack */
  ReliabilityHelper reliabilityHelper;
  // Configure the power, temperature and reliability model for each node
  for (uint8_t i = 0; i < nDevices; ++i)
  {
    Ptr<Node> nodei = sensorDevices.Get(i);
    reliabilityHelper.SetDeviceType("Arduino");
    reliabilityHelper.SetPowerModel("ns3::StatePowerModel",
      "IdlePowerW", DoubleValue (P0 + SensorFlag[i] * Es / packetInterval), 
      "TxPowerW", DoubleValue(P0 + SensorFlag[i] * Es / packetInterval + 1*ptxo[i]),
      "TxDurationS", DoubleValue(originalPacketSize / bw),
      "RxPowerW", DoubleValue(P0 + SensorFlag[i] * Es / packetInterval + Prx),
      "RxDurationS", DoubleValue(originalPacketSize / bw));
    reliabilityHelper.SetPerformanceModel("ns3::PerformanceSimpleModel");
    reliabilityHelper.SetTemperatureModel("ns3::TemperatureSimpleModel",
      "SimSpeed", IntegerValue(simSpeed),
      "A", DoubleValue(1.0),
      "B", DoubleValue(0.0),
      "C", DoubleValue(100.0),
      "D", DoubleValue(0.0));
    reliabilityHelper.SetAmbientTemperature(tempMatrix[i]);
    reliabilityHelper.SetReliabilityModel("ns3::ReliabilityTDDBModel",
      "SimSpeed", IntegerValue(simSpeed));
    // reliabilityHelper.SetApplication("AdaBoost",dataSize,packetSize);
    reliabilityHelper.Install(nodei);
  }
  
  /* cpu energy model */
  CpuEnergyModelHelper cpuEnergyHelper;
  cpuEnergyHelper.SetSimSpeed(simSpeed);
  DeviceEnergyModelContainer deviceModels = cpuEnergyHelper.Install(devicesNet, energySources);


  /* energy harvester */
  BasicEnergyHarvesterHelper basicHarvesterHelper;
  // configure energy harvester
  basicHarvesterHelper.Set ("SimSpeed", IntegerValue(simSpeed));
  basicHarvesterHelper.Set ("PeriodicHarvestedPowerUpdateInterval", TimeValue (Seconds (harvestingUpdateInterval)));
  basicHarvesterHelper.Set ("HarvestablePower", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=0.1]"));
  // basicHarvesterHelper.SetRechargeTrace (rcMatrix[i]);
  // install harvester on all energy sources
  EnergyHarvesterContainer harvesters = basicHarvesterHelper.Install (energySources, rcMatrix);
  // std::cout<<"Energy harvesters installed.\n";
  /***************************************************************************/

  /** connect trace sources **/
  /***************************************************************************/
  // // energy source
  Ptr<BasicEnergySource> basicSourcePtr = DynamicCast<BasicEnergySource> (energySources.Get (6));
  // basicSourcePtr->TraceConnectWithoutContext ("RemainingEnergy", MakeCallback (&RemainingEnergy));
  // device energy model
  Ptr<DeviceEnergyModel> basicRadioModelPtr =
    basicSourcePtr->FindDeviceEnergyModels ("ns3::CpuEnergyModel").Get (0);
  NS_ASSERT (basicRadioModelPtr != NULL);
  //  basicRadioModelPtr->TraceConnectWithoutContext ("TotalEnergyConsumption", MakeCallback (&AveragePower));
  //  basicRadioModelPtr->TraceConnectWithoutContext ("TotalEnergyConsumption", MakeCallback (&TotalEnergy));
  // energy harvester
  Ptr<BasicEnergyHarvester> basicHarvesterPtr = DynamicCast<BasicEnergyHarvester> (harvesters.Get (6));
  // basicHarvesterPtr->TraceConnectWithoutContext ("HarvestedPower", MakeCallback (&HarvestedPower));
  //  basicHarvesterPtr->TraceConnectWithoutContext ("TotalEnergyHarvested", MakeCallback (&TotalEnergyHarvested));
  // /***************************************************************************/

  /////////////
  // Tracing //
  /////////////
  if (tracing == true)
  {
    AsciiTraceHelper ascii;
    wifiPhy.EnableAsciiAll (ascii.CreateFileStream ("reliability-example.tr"));
    wifiPhy.EnablePcap ("reliability-example", allNodes);

    // Trace routing tables
    Ptr<OutputStreamWrapper> routingStream = Create<OutputStreamWrapper> ("reliability-example.routes", std::ios::out);
    staticRouting.PrintRoutingTableAllEvery (Seconds (2), routingStream);
    Ptr<OutputStreamWrapper> neighborStream = Create<OutputStreamWrapper> ("reliability-example.neighbors", std::ios::out);
    staticRouting.PrintNeighborCacheAllEvery (Seconds (2), neighborStream);

  }

  // Configure the power, temperature and reliability model for each node
  // for (uint8_t i = 0; i < nDevices; ++i)
  // {
  //   Ptr<Node> nodei = sensorDevices.Get(i);
  //   PrintInfo(nodei);
  // }
  // double energyPeriod = 1; // print period (hours)
  // for (uint8_t i = 0; i < nDevices; ++i)
  // {
  //   Ptr<EnergySource> sourcei = energySources.Get(i);
  //   PrintEnergy(sourcei,energyPeriod*3600/simSpeed); // hours to simulation time in seconds
  // }

  // double reliabilityPeriod = 24*7; // print period (hours)
  // for (uint8_t i = 0; i < nDevices; ++i)
  // {
  //   Ptr<Node> nodei = sensorDevices.Get(i);
  //   PrintReliability(nodei,reliabilityPeriod*3600/simSpeed); // hours to simulation time in seconds
  // }

  // double energyPeriod = 1; // print period (hours)
  // PrintEnergy(energyResults,energySources,energyPeriod*3600/simSpeed); // hours to simulation time in seconds

  double energyPeriod = 1; // print period (hours)
  PrintEnergy(engFile,energySources,energyPeriod*3600/simSpeed); // hours to simulation time in seconds


  // int schedulingIdx = 1;
  // while (schedulingIdx*energyPeriod<simulationSeconds)
  // {
  //   Simulator::Schedule (Seconds (schedulingIdx*energyPeriod), &PrintEnergy,energyResults,energySources,energyPeriod*3600/simSpeed);
  //   schedulingIdx++;
  // }

  // double reliabilityPeriod = 24*7; // print period (hours)
  // PrintReliability(reliabilityResults,sensorDevices,reliabilityPeriod*3600/simSpeed); // hours to simulation time in seconds

  double reliabilityPeriod = 30*24/5; // print period (hours)
  PrintReliability(relFile,sensorDevices,reliabilityPeriod*3600/simSpeed); // hours to simulation time in seconds
 
  // schedulingIdx = 1;
  // while (schedulingIdx*reliabilityPeriod<simulationSeconds)
  // {
  //   Simulator::Schedule (Seconds (schedulingIdx*reliabilityPeriod), &PrintReliability,reliabilityResults,sensorDevices,reliabilityPeriod*3600/simSpeed);
  //   schedulingIdx++;
  // }  


  Simulator::Schedule (Seconds (simulationSeconds-1), &PrintAveragePower,pwrFile, energySources,simSpeed);

  // PrintTime();

  /** Flow Monitor **/
  /***************************************************************************/
  Ptr<FlowMonitor> flowMonitor;
  FlowMonitorHelper flowHelper;
  if (flow_monitor)
    {
      flowMonitor = flowHelper.InstallAll ();
    }

  /** simulation setup **/
  Simulator::Stop (Seconds (simulationSeconds));
  //Simulator::Stop (Seconds (100));
  Simulator::Run ();

  if (flow_monitor)
    {
      flowMonitor->SerializeToXmlFile ("flowmonitor.xml", true, true);
    }

  // Print per flow statistics
  flowMonitor->CheckForLostPackets ();
  Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowHelper.GetClassifier ());
  FlowMonitor::FlowStatsContainer stats = flowMonitor->GetFlowStats ();
  Histogram delayHistogram;
  for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin (); i != stats.end (); ++i)
    {
      Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (i->first);
      

          std::cout << "Flow " << i->first << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
          std::cout << "  Tx Packets: " << i->second.txPackets << "\n";
          std::cout << "  Rx Packets: " << i->second.rxPackets << "\n";
          std::cout << "  DelaySum: " << i->second.delaySum.GetSeconds() << "\n";
          std::cout << "  Mean Delay: " << (i->second.delaySum.GetSeconds())/(i->second.rxPackets) << "\n";
          delayHistogram = i->second.delayHistogram;       
          std::cout << "  Number of Hops: " << i->second.timesForwarded << "\n";   
          std::cout << "  Throughput: " << i->second.rxBytes / (i->second.timeLastRxPacket.GetSeconds()-i->second.timeFirstTxPacket.GetSeconds()) << " Bytes/s\n";
    }
  /***************************************************************************/


  std::ofstream sohResults;
  sohResults.open(sohFile);
  std::ofstream mttfResults;
  mttfResults.open(mttfFile);

  
// SoH and MTTF calculation
  std::vector<double> soh_vec;
  std::vector<double> mttf_vec;
  for (uint8_t i = 0; i < nDevices; ++i)
  {
    Ptr<Node> nodei = sensorDevices.Get(i);
    soh_vec.push_back(compute_soh(nodei->GetObject<TemperatureModel>()->GetTemperatureVector()));
    mttf_vec.push_back(compute_mttf(nodei->GetObject<TemperatureModel>()->GetTemperatureVector()));
  }
  for (uint8_t i = 0; i < nDevices; ++i)
  {
    std::cout<< "SoH " << unsigned(i) << ": "<<soh_vec[i]<<"\n";
    sohResults<< soh_vec[i] <<"\n";
  }
  for (uint8_t i = 0; i < nDevices; ++i)
  {
    std::cout<< "MTTF " << unsigned(i) << ": "<<mttf_vec[i]<<"\n";
    mttfResults<< mttf_vec[i] <<"\n";
  } 
  sohResults.close ();
  mttfResults.close ();

  // energyResults.close ();
  // reliabilityResults.close ();

  Simulator::Destroy ();

  return 0;
}

// split implementation for reading from external text files
std::vector < std::string > split(std::string const & str, const char delim)
{
    std::vector < std::string > result;

    std::stringstream ss(str);
    std::string s;

    while (std::getline(ss, s, delim)) {
        result.push_back(s);
    }

    return result;
}


// SoH computation
double compute_soh(std::vector<double> Tc)
{
 int N_bin = 20;
 double histogrammdistance = (ceil(*max_element(Tc.begin(), Tc.end()))-floor(*min_element(Tc.begin(), Tc.end())))/N_bin;
 std::vector<double> edges;
 double min_temp = *min_element(Tc.begin(), Tc.end());

 for (int i = 0 ;i <=N_bin; i++)
 {
     edges.push_back(floor(min_temp) +  histogrammdistance * i);
 }
 std::vector<double> counts(N_bin);

 for (uint32_t i = 0 ;i <unsigned(Tc.size()); i++)
 {
  for (int j = 0 ;j <=N_bin-1; j++)
  {
      if((Tc[i]>edges[j])&&(Tc[i]<edges[j+1]))
      {
        counts[j] = counts[j] + 1.0/Tc.size();
      }
  }
 }
 double Tref = 25 + 273.15;
 double Tsec  = 5 * 365 * 24 *3600;
 double kt = 4.14e-10;
 double kT = 6.93e-2;
 double SoH = 0;
 for (int i = 0 ;i <=N_bin-1; i++)
  {
      SoH = SoH + counts[i]*exp(-kt * Tsec * exp(kT * Tref * (1 - (Tref / ((edges[i]+edges[i+1])/2)))));
  }
  return SoH;
 
}


// MTTF computation
double compute_mttf(std::vector<double> Tc)
{
 int N_bin = 20;
 double histogrammdistance = (ceil(*max_element(Tc.begin(), Tc.end()))-floor(*min_element(Tc.begin(), Tc.end())))/N_bin;
 std::vector<double> edges;
 double min_temp = *min_element(Tc.begin(), Tc.end());

 for (int i = 0 ;i <=N_bin; i++)
 {
     edges.push_back(floor(min_temp) +  histogrammdistance * i);
 }
 std::vector<double> counts(N_bin);

 for (uint32_t i = 0 ;i <unsigned(Tc.size()); i++)
 {
  // std::cout << edges[i] << ' ';
  for (int j = 0 ;j <=N_bin-1; j++)
  {
      if((Tc[i]>edges[j])&&(Tc[i]<edges[j+1]))
      {
        counts[j] = counts[j] + 1.0/Tc.size();
      }
  }
 }
 double Tref = 25 + 273.15;
 double c  = 10 * 1.1949e3 / (6.022 * 1.38);
 double mttf = 0;
 for (int i = 0 ;i <=N_bin-1; i++)
  {
      mttf = mttf + counts[i]*exp(c / ((edges[i]+edges[i+1])/2)) / exp(c / Tref);
  }
  
  return mttf;
}