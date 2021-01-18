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

#ifndef POWER_LINEARMODEL_H
#define POWER_LINEARMODEL_H

#include "ns3/nstime.h"
#include "ns3/event-id.h"
#include "ns3/traced-value.h"
#include "ns3/temperature-model.h"
#include <ns3/performance-model.h>
#include "ns3/power-model.h"
#include "ns3/wifi-phy.h"


namespace ns3 {

class PowerDistModel : public PowerModel
{
public:

  static TypeId GetTypeId (void);
  PowerDistModel ();

  virtual ~PowerDistModel ();



  /**
   * \param  Pointer to temperature object attached to the device.
   *
   * Registers the Temperature Model to Power Model.
   */
  virtual void RegisterTemperatureModel (Ptr<TemperatureModel> temperatureModel);
  /**
   * \param  Pointer to performance object attached to the device.
   *
   * Registers the Performance Model to Power Model.
   */
  virtual void RegisterPerformanceModel (Ptr<PerformanceModel> performanceModel);
  /**
   * \brief Sets pointer to EnergySouce installed on node.
   *
   * \param source Pointer to EnergySource installed on node.
   *
   * Implements DeviceEnergyModel::SetEnergySource.
   */
  //virtual void SetEnergySource (Ptr<EnergySource> source);

  // Setter & getters.
  virtual double GetIdlePowerW (void) const;
  virtual void SetIdlePowerW (double IdlePowerW);
  virtual double GetTxPowerW (void) const;
  virtual void SetTxPowerW (double TxPowerW);
  virtual double GetTxDurationS (void) const;
  virtual void SetTxDurationS (double TxDurationS);
  virtual double GetRxPowerW (void) const;
  virtual void SetRxDurationS (double RxDurationS);
  virtual double GetIdlePowerW (void) const;
  virtual void SetIdlePowerW (double IdlePowerW);
  virtual int GetState (void) const;
  virtual void SetState (int state);
  
  /**
   * \returns Current power.
   */
  virtual double GetPower (void) const;

  /**
   * \returns Total energy to be consumed.
   */
  virtual double GetEnergy (void) const;

  /**
   * \brief Updates the power.
   *
   */
  virtual void UpdatePower ();


private:
  virtual void DoDispose (void);

private:

  Ptr<TemperatureModel> m_temperatureModel;
  Ptr<PerformanceModel> m_performanceModel;

  // Member variables for current draw in different radio modes.
  WifiPhy::State m_currentState;
  double m_idlePowerW;
  double m_txPowerW;
  double m_txDurationS;
  double m_rxPowerW;
  double m_rxDurationsS;
  // This variable keeps track of the total energy consumed by this model.
  TracedValue<double> m_cpupower;
  // State variables.
  EventId m_powerUpdateEvent;            // energy update event
  Time m_lastUpdateTime;          // time stamp of previous energy update
  Time m_powerUpdateInterval;            // energy update interval

};

} // namespace ns3

#endif /* POWER_DIST_MODEL_H */
