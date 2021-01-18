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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("ReliabilityExample");

static inline std::string
PrintReceivedPacket (Address& from)
{
  InetSocketAddress iaddr = InetSocketAddress::ConvertFrom (from);

  std::ostringstream oss;
  oss << "--\nReceived one packet! Socket: " << iaddr.GetIpv4 ()
      << " port: " << iaddr.GetPort ()
      << " at time = " << Simulator::Now ().GetSeconds ()
      << "\n--";

  return oss.str ();
}

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
          NS_LOG_UNCOND (PrintReceivedPacket (from));
        }
    }
}

/**
 * \param socket Pointer to socket.
 * \param pktSize Packet size.
 * \param n Pointer to node.
 * \param pktCount Number of packets to generate.
 * \param pktInterval Packet sending interval.
 *
 * Traffic generator.
 */
static void
GenerateTraffic (Ptr<Socket> socket, uint32_t pktSize, Ptr<Node> n,
                 uint32_t pktCount, Time pktInterval)
{
  if (pktCount > 0)
    {
      socket->Send (Create<Packet> (pktSize));
      Simulator::Schedule (pktInterval, &GenerateTraffic, socket, pktSize, n,
                           pktCount - 1, pktInterval);
    }
  else
    {
      socket->Close ();
    }
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

  std::cout<<"At time "<< Simulator::Now().GetSeconds()<<", NodeId = "<<node->GetId();
  std::cout << " CPU Power = " << node->GetObject<PowerModel>()->GetPower();
  std::cout << " Performance = " << node->GetObject<PerformanceModel>()->GetThroughput();
  std::cout << " Temperature = " << node->GetObject<TemperatureModel>()->GetTemperature()<<std::endl;
  std::cout << " Reliability = " << node->GetObject<ReliabilityModel>()->GetReliability()<<std::endl;

  if (!Simulator::IsFinished ())
  {
    Simulator::Schedule (Seconds (0.5),&PrintInfo,node);
  }
}

std::string srlocFile = "dev_loc.txt"; // Sensor location file
std::string gwlocFile = "gw_loc.txt"; // Gateway location file
std::string flFile = "fl_mat.txt"; // Flow quantity file
std::string tempFile = "temp.txt"; // Temperature traces at each sensor location


int
main (int argc, char *argv[])
{
  /*
  LogComponentEnable ("CpuEnergyModel", LOG_LEVEL_DEBUG);
  LogComponentEnable ("PowerModel", LOG_LEVEL_DEBUG);
  LogComponentEnable ("TemperatureSimpleModel", LOG_LEVEL_DEBUG);
  LogComponentEnable ("ReliabilityTDDBModel", LOG_LEVEL_DEBUG);
   */
  LogComponentEnable ("AppPowerModel", LOG_LEVEL_DEBUG);
  
 
  int nDevices = 0;
  int nGateways = 0;
  std::string phyMode ("DsssRate1Mbps");
  double Prss = -80;            // dBm
  uint32_t PpacketSize = 100;   // bytes
  bool verbose = false;
  uint32_t dataSize = 10000;    // bytes for reliability helper
  // simulation parameters
  uint32_t numPackets = 10000;  // number of packets to send
  double interval = 10;          // seconds
  double startTime = 0.0;       // seconds
  /*
   * This is a magic number used to set the transmit power, based on other
   * configuration.
   */
  double offset = 81;

  CommandLine cmd;
  cmd.AddValue ("phyMode", "Wifi Phy mode", phyMode);
  cmd.AddValue ("Prss", "Intended primary RSS (dBm)", Prss);
  cmd.AddValue ("PpacketSize", "size of application packet sent", PpacketSize);
  cmd.AddValue ("numPackets", "Total number of packets to send", numPackets);
  cmd.AddValue ("startTime", "Simulation start time", startTime);
  cmd.AddValue ("distanceToRx", "X-Axis distance between nodes", distanceToRx);
  cmd.Parse (argc, argv);

  // Convert to time object
  Time interPacketInterval = Seconds (interval);

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
  std::ifstream EdLocationFile(srlocFile);
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
    NS_LOG_ERROR ("Unable to open file " << srlocFile);
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
  std::ifstream GwLocationFile(gwlocFile);
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
    NS_LOG_ERROR ("Unable to open file " << gwlocFile);
    return -1;
  }

  gateways.Create (nGateways);
  mobilityGw.Install (gateways);

  NodeContainer allNodes = NodeContainer(sensorDevices, gateways);


  ////////////////
  // Wifi PHY   //
  ////////////////

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

  // Add an upper mac and disable rate control
  WifiMacHelper wifiMac;
  wifi.SetStandard (WIFI_PHY_STANDARD_80211b);
  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
                                "DataMode",StringValue (phyMode),
                                "ControlMode",StringValue (phyMode));
  // Set it to adhoc mode
  wifiMac.SetType ("ns3::AdhocWifiMac");

  /** install PHY + MAC **/
  NetDeviceContainer devicesNet = wifi.Install (wifiPhy, wifiMac, sensorDevices);
  NetDeviceContainer gatewaysNet = wifi.Install (wifiPhy, wifiMac, gateways);
  NetDeviceContainer allNodesNet = NetDeviceContainer(devicesNet, gatewaysNet);


  /** Energy Model **/
  /***************************************************************************/
  /* energy source */
  //BasicEnergySourceHelper basicSourceHelper;
  // configure energy source
  //basicSourceHelper.Set ("BasicEnergySourceInitialEnergyJ", DoubleValue (100000));
  // install source
  //EnergySourceContainer source1 = basicSourceHelper.Install (node1);
  /* reliability stack */
  //ReliabilityHelper reliabilityHelper;
  //reliabilityHelper.SetDeviceType("RaspberryPi");
  //reliabilityHelper.SetPowerModel("ns3::AppPowerModel");
  //reliabilityHelper.SetPerformanceModel("ns3::PerformanceSimpleModel");
  //reliabilityHelper.SetTemperatureModel("ns3::TemperatureSimpleModel");
  //reliabilityHelper.SetReliabilityModel("ns3::ReliabilityTDDBModel");
  //reliabilityHelper.SetApplication("AdaBoost",dataSize,PpacketSize);
  //reliabilityHelper.Install(node1);
  /* cpu energy model */
  // CpuEnergyModelHelper cpuEnergyHelper;
  // DeviceEnergyModelContainer deviceModels = cpuEnergyHelper.Install(device1, source1);
  /***************************************************************************/


  /** Internet stack **/
  InternetStackHelper internet;
  internet.Install (allNodes);

  Ipv4AddressHelper ipv4;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer allNodesInterface = ipv4.Assign (allNodesNet);

  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  for (int i = 0; i < nGateways; ++i)
  {
    Ptr<Socket> recvSink = Socket::CreateSocket (gateways.Get(i), tid);  // only one sink
    InetSocketAddress local = InetSocketAddress (Ipv4Address::GetAny (), 80);
    recvSink->Bind (local);
    recvSink->SetRecvCallback (MakeCallback (&ReceivePacket));
  }

  Ptr<Socket> source = Socket::CreateSocket (node0, tid);    // node 0, sender
  InetSocketAddress remote = InetSocketAddress (Ipv4Address::GetBroadcast (), 80);
  source->SetAllowBroadcast (true);
  source->Connect (remote);


  PrintInfo (node1);


  /** simulation setup **/
  // start traffic
  Simulator::Schedule (Seconds (startTime), &GenerateTraffic, source, PpacketSize,
                       networkNodes.Get (0), numPackets, interPacketInterval);

  Simulator::Stop (Seconds (100.0));
  Simulator::Run ();
  Simulator::Destroy ();

  return 0;
}
