/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
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

#include "ns3/log.h"
#include "ns3/traced-value.h"
#include "ns3/double.h"
#include "ns3/simulator.h"
#include "ns3/trace-source-accessor.h"
#include "ns3/pointer.h"
#include "ns3/power-model.h"
#include "ns3/power-distmodel.h"
#include <ns3/performance-model.h>
#include <iterator>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream> //istringstream


NS_LOG_COMPONENT_DEFINE ("PowerDistModel");

namespace ns3 {

NS_OBJECT_ENSURE_REGISTERED (PowerDistModel);

TypeId
PowerDistModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::PowerDistModel")
    .SetParent<PowerModel> ()
    .SetGroupName ("Power")
    .AddConstructor<PowerDistModel> ()
    .AddAttribute ("TxPowerW",
                   "Transmission Power of WiFi Chip in W",
                   DoubleValue (0.22),    // default
                   MakeDoubleAccessor (&PowerDistModel::SetTxPowerW,
                                       &PowerDistModel::GetTxPowerW),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("TxDurationS",
                   "Transmission duration of WiFi Chip in S",
                   DoubleValue (0),    // default
                   MakeDoubleAccessor (&PowerDistModel::SetTxDurationS,
                                       &PowerDistModel::GetTxDurationS),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("RxPowerW",
                   "Receiption Power of WiFi Chip in W",
                   DoubleValue (0.1),    // default
                   MakeDoubleAccessor (&PowerDistModel::SetRxPowerW,
                                       &PowerDistModel::GetRxPowerW),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("RxDurationS",
                   "Receiption duration of WiFi Chip in S",
                   DoubleValue (0),    // default
                   MakeDoubleAccessor (&PowerDistModel::SetRxDurationS,
                                       &PowerDistModel::GetRxDurationS),
                   MakeDoubleChecker<double> ())
    .AddAttribute ("IdlePowerW",
                   "Idle Power Consumption of Cpu",
                   DoubleValue (2.8),    // default
                   MakeDoubleAccessor (&PowerDistModel::SetIdlePowerW,
                                       &PowerDistModel::GetIdlePowerW),
                   MakeDoubleChecker<double> ())
    .AddTraceSource ("CpuPower",
                     "CPU power consumption of the device.",
                     MakeTraceSourceAccessor (&PowerDistModel::m_cpupower),
                     "ns3::TracedValueCallback::Double")
  ; 
  return tid;
}

PowerDistModel::PowerDistModel ()
{
  NS_LOG_FUNCTION (this);
  m_lastUpdateTime = Seconds (0.0);
  m_powerUpdateInterval = Seconds (0.045);
  m_temperatureModel = NULL;      // TemperatureModel
  m_performanceModel = NULL;
  m_currentState = 0;
}

PowerDistModel::~PowerDistModel ()
{
  NS_LOG_FUNCTION (this);
}


void
PowerDistModel::RegisterTemperatureModel (Ptr<TemperatureModel> temperatureModel)
{
  m_temperatureModel = temperatureModel;
}

void
PowerDistModel::RegisterPerformanceModel (Ptr<PerformanceModel> performanceModel)
{
  m_performanceModel = performanceModel;
}

double
PowerDistModel::GetPower (void) const
{
  NS_LOG_FUNCTION (this);
  return m_cpupower;
}

double
PowerDistModel::GetDuration (void) const
{
  NS_LOG_FUNCTION (this);
  return m_duration;
}

void
PowerDistModel::SetIdlePowerW (double IdlePowerW)
{
  NS_LOG_FUNCTION (this << IdlePowerW);
  m_idlePowerW = IdlePowerW;
}

double
PowerDistModel::GetIdlePowerW (void) const
{
  NS_LOG_FUNCTION (this);
  return m_idlePowerW;
}

void
PowerDistModel::SetTxPowerW (double TxPowerW)
{
  NS_LOG_FUNCTION (this << TxPowerW);
  m_txPowerW = TxPowerW;
}

double
PowerDistModel::GetTxPowerW (void) const
{
  NS_LOG_FUNCTION (this);
  return m_txPowerW;
}

void
PowerDistModel::SetTxDurationS (double TxDurationS)
{
  NS_LOG_FUNCTION (this << TxDurationS);
  m_txDurationS = TxDurationS;
}

double
PowerDistModel::GetTxDurationS (void) const
{
  NS_LOG_FUNCTION (this);
  return m_txDurationS;
}

void
PowerDistModel::SetRxPowerW (double RxPowerW)
{
  NS_LOG_FUNCTION (this << RxPowerW);
  m_rxPowerW = RxPowerW;
}

double
PowerDistModel::GetRxPowerW (void) const
{
  NS_LOG_FUNCTION (this);
  return m_rxPowerW;
}

void
PowerDistModel::SetRxDurationS (double RxDurationS)
{
  NS_LOG_FUNCTION (this << RxDurationS);
  m_rxDurationS = RxDurationS;
}

double
PowerDistModel::GetRxDurationS (void) const
{
  NS_LOG_FUNCTION (this);
  return m_rxDurationS;
}

void
PowerDistModel::SetState (int state)
{
  NS_LOG_FUNCTION (this << state);
  m_currentState = state;
}

int
PowerDistModel::GetState (void) const
{
  NS_LOG_FUNCTION (this);
  return m_currentState;
}

void
PowerDistModel::UpdatePower ()
{
  NS_LOG_FUNCTION (this);
  NS_LOG_DEBUG ("PowerDistModel:Updating power" << " at time = " << Simulator::Now ());
  if (Simulator::IsFinished ())
  {
    return;
  }
  m_powerUpdateEvent.Cancel ();

  switch(m_currentState)
  {
    case WifiPhy::IDLE:
      m_cpupower = m_idlePowerW;
      m_duration = 0.0; // does not matter, will not be used
      break;
    case WifiPhy::TX:
      m_cpupower = m_txPowerW;
      m_duration = m_txDurationS;
      break;
    case WifiPhy::RX:
      m_cpupower = m_rxPowerW;
      m_duration = m_rxDurationS;
      break;
    default:
      m_cpupower = m_idlePowerW;
      m_duration = 0.0; // does not matter, will not be used
  }

  // update last update time stamp
  m_lastUpdateTime = Simulator::Now ();

  m_temperatureModel->UpdateTemperature (m_cpupower);
  m_powerUpdateEvent = Simulator::Schedule (m_powerUpdateInterval,&PowerDistModel::UpdatePower,this);
}

void
PowerDistModel::DoDispose (void)
{
  NS_LOG_FUNCTION (this);
  m_temperatureModel = NULL;      // TemperatureModel
  m_performanceModel = NULL;      // PerformanceModel

}



} // namespace ns3
