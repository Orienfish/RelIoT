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
#include "ns3/internet-module.h"
#include "ns3/olsr-helper.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("ReliabilityExample");

std::vector < std::string > split(std::string const & str, const char delim);

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

/// Trace function for total energy consumption at node.
void
TotalEnergy (double oldValue, double totalEnergy)
{
  NS_LOG_UNCOND (Simulator::Now ().GetSeconds ()
                 << "s Total energy consumed by radio = " << totalEnergy << "J");
}

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
    Simulator::Schedule (Seconds (0.5),&PrintInfo,node);
  }
}

std::string srFile = "dev.txt";     // Sensor location file
std::string gwFile = "gw.txt";      // Gateway location file
std::string flFile = "fl.txt";      // File of transmission distance at each node
std::string distFile = "dist.txt";  // File of transmission distance at each node
std::string tempFile = "temp.txt";  // Temperature traces at each sensor location


int
main (int argc, char *argv[])
{
  // LogComponentEnable ("CpuEnergyModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("PowerModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("TemperatureSimpleModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("ReliabilityTDDBModel", LOG_LEVEL_DEBUG);
  // LogComponentEnable ("AppPowerModel", LOG_LEVEL_DEBUG);

  uint32_t nDevices = 0;
  uint32_t nGateways = 0;
  std::string phyMode ("DsssRate1Mbps");
  double Prss = -80;            // dBm
  uint32_t packetSize = 100;   // bytes
  double bw = 2000.0;             // B/s
  // uint32_t dataSize = 10000;    // bytes for reliability helper

  // Energy parameters
  double Es = 0.04;             // J, energy for one sensing sample
  double P0 = 0.01;             // W, ambient power dissipation
  double Pto = 0.22;            // W, ambient wifi transmission power
  double alpha = 3.5;           // Path exponent
  double beta = 0.0000001;      // Path linear coefficient
  double Prx = 0.1;             // W, wifi receiption power
  double BattJ = 23760;         // J, battery capacity

  // Simulation parameters
  double packetInterval = 5;    // seconds
  double startTime = 0.0;       // seconds
  double simulationTime = 1.0;  // years
  uint32_t port = 5000;         // port number that receiver is listening on
  bool verbose = false;
  bool tracing = false;

  /*
   * This is a magic number used to set the transmit power, based on other
   * configuration.
   */
  double offset = 81;

  CommandLine cmd;
  cmd.AddValue ("phyMode", "Wifi Phy mode", phyMode);
  cmd.AddValue ("Prss", "Intended primary RSS (dBm)", Prss);
  cmd.AddValue ("packetSize", "size of application packet sent", packetSize);
  cmd.AddValue ("startTime", "Simulation start time", startTime);
  cmd.Parse (argc, argv);

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


  //////////////////////////
  // Read distance matrix //
  //////////////////////////
  std::ifstream DistFile(distFile);
  std::vector<double> dist; // flow matrix
  if (DistFile.is_open())
  {
    NS_LOG_DEBUG ("Read from existing distance file.");
    std::string line;
    while (std::getline(FLFile, line)) {
        if (line.size() > 0) {
            std::vector < std::string > data = split(line, ' ');
            dist.push_back(atof(data.at(0).c_str()));
        }
    }
  }
  else
  {
    NS_LOG_ERROR ("Unable to open file " << distFile);
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
  internet.SetRoutingHelper (staticRouting); // has effect on the next Install ()
  for (int i = 0; i < nDevices; ++i)
  {
    for (int j = 0; j < nDevices + nGateways; ++j)
    {
      if (fij[i][j] < 0.1) // If no flow is assigned, skip this connection
        continue;

      Ptr<Node>srNode = sensorDevices.Get(i);
      Ptr<Ipv4> ipv4Node = srNode->GetObject<Ipv4> ();
      Ptr<Ipv4StaticRouting> staticRoutingNode = staticRouting.GetStaticRouting (ipv4Node);
      Ipv4Address sinkAddr = allNodesInterface.GetAddress(nDevices+nGateways-1);
      Ipv4Address nextHopAddr = allNodesInterface.GetAddress(j);
      staticRoutingNode->AddHostRouteTo (sinkAddr, nextHopAddr, 1);
      sinkAddr.Print(std::cout);
      nextHopAddr.Print(std::cout);
    }
  }

  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  // Configure receiving node - the only sink
  Ipv4Address sinkAddr = allNodesInterface.GetAddress(nDevices+nGateways-1);
  Ptr<Socket> recvSocket = Socket::CreateSocket (allNodes.Get(nDevices+nGateways-1), tid);
  InetSocketAddress local = InetSocketAddress (sinkAddr, port);
  recvSocket->Bind (local);
  recvSocket->SetRecvCallback (MakeCallback (&ReceivePacket));

  // Configure source node
  for (int i = 0; i < nDevices; ++i)
  {
    if (SensorFlag[i] > 0)
    {
      Ptr<Socket> source = Socket::CreateSocket (sensorDevices.Get(i), tid);
      InetSocketAddress remote = InetSocketAddress (sinkAddr, port);
      source->Connect (remote);

      // Schedule SendPacket
      Simulator::ScheduleWithContext(source->GetNode()->GetId(),
                                     Seconds (startTime), &SendPacket,
                                     source, packetSize,
                                     Seconds(packetInterval));
    }
  }


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
  for (int i = 0; i < nDevices; ++i)
  {
    Ptr<Node> nodei = sensorDevices.Get(i);
    reliabilityHelper.SetDeviceType("Arduino");
    reliabilityHelper.SetPowerModel("ns3::PowerDistModel",
      "IdlePower", DoubleValue (P0 + SensorFlag[i] * Es / packetInterval), 
      "TxPowerW", DoubleValue(Pto + beta * pow(dist[i], alpha)),
      "TxDurationS", DoubleValue(packetSize / bw),
      "RxPowerW", DoubleValue(Prx),
      "RxDurationS", DoubleValue(packetSize / bw));
    reliabilityHelper.SetPerformanceModel("ns3::PerformanceSimpleModel");
    reliabilityHelper.SetTemperatureModel("ns3::TemperatureSimpleModel");
    reliabilityHelper.SetReliabilityModel("ns3::ReliabilityTDDBModel");
    // reliabilityHelper.SetApplication("AdaBoost",dataSize,packetSize);
    reliabilityHelper.Install(nodei);
  }
  
  /* cpu energy model */
  CpuEnergyModelHelper cpuEnergyHelper;
  DeviceEnergyModelContainer deviceModels = cpuEnergyHelper.Install(devicesNet, energySources);

  /***************************************************************************/
  if (tracing == true)
  {
    AsciiTraceHelper ascii;
    wifiPhy.EnableAsciiAll (ascii.CreateFileStream ("reliability-example.tr"));
    wifiPhy.EnablePcap ("reliability-example", allNodes);
  }

  /** simulation setup **/
  // Simulator::Stop (Years (simulationTime));
  Simulator::Stop (Seconds (100));
  Simulator::Run ();
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
