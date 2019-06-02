/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2017 Vishwesh Rege
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

#include "reliability-helper.h"
#include <ns3/log.h>
#include "ns3/names.h"
#include <vector>
namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("ReliabilityHelper");

ReliabilityHelper::ReliabilityHelper (void)
{
  m_power.SetTypeId ("ns3::PowerLinearModel");
  m_performance.SetTypeId ("ns3::PerformanceSimpleModel");
  m_temperature.SetTypeId ("ns3::TemperatureModel");
  m_appname = "LinearRegression";
  m_datasize = 1000;
}

ReliabilityHelper::~ReliabilityHelper (void)
{
}

 
void 
ReliabilityHelper::SetPowerModel (std::string type,
                                  std::string n1, const AttributeValue &v1,
                                  std::string n2, const AttributeValue &v2,
                                  std::string n3, const AttributeValue &v3,
                                  std::string n4, const AttributeValue &v4,
                                  std::string n5, const AttributeValue &v5,
                                  std::string n6, const AttributeValue &v6,
                                  std::string n7, const AttributeValue &v7,
                                  std::string n8, const AttributeValue &v8,
                                  std::string n9, const AttributeValue &v9)
{
  m_power.SetTypeId (type);
  m_power.Set (n1, v1);
  m_power.Set (n2, v2);
  m_power.Set (n3, v3);
  m_power.Set (n4, v4);
  m_power.Set (n5, v5);
  m_power.Set (n6, v6);
  m_power.Set (n7, v7);
  m_power.Set (n8, v8);
  m_power.Set (n9, v9);
}

void 
ReliabilityHelper::SetPerformanceModel (std::string type,
                                  std::string n1, const AttributeValue &v1,
                                  std::string n2, const AttributeValue &v2,
                                  std::string n3, const AttributeValue &v3,
                                  std::string n4, const AttributeValue &v4,
                                  std::string n5, const AttributeValue &v5,
                                  std::string n6, const AttributeValue &v6,
                                  std::string n7, const AttributeValue &v7,
                                  std::string n8, const AttributeValue &v8,
                                  std::string n9, const AttributeValue &v9)
{
  m_performance.SetTypeId (type);
  m_performance.Set (n1, v1);
  m_performance.Set (n2, v2);
  m_performance.Set (n3, v3);
  m_performance.Set (n4, v4);
  m_performance.Set (n5, v5);
  m_performance.Set (n6, v6);
  m_performance.Set (n7, v7);
  m_performance.Set (n8, v8);
  m_performance.Set (n9, v9);
}

void 
ReliabilityHelper::SetTemperatureModel (std::string type,
                                  std::string n1, const AttributeValue &v1,
                                  std::string n2, const AttributeValue &v2,
                                  std::string n3, const AttributeValue &v3,
                                  std::string n4, const AttributeValue &v4,
                                  std::string n5, const AttributeValue &v5,
                                  std::string n6, const AttributeValue &v6,
                                  std::string n7, const AttributeValue &v7,
                                  std::string n8, const AttributeValue &v8,
                                  std::string n9, const AttributeValue &v9)
{
  m_temperature.SetTypeId (type);
  m_temperature.Set (n1, v1);
  m_temperature.Set (n2, v2);
  m_temperature.Set (n3, v3);
  m_temperature.Set (n4, v4);
  m_temperature.Set (n5, v5);
  m_temperature.Set (n6, v6);
  m_temperature.Set (n7, v7);
  m_temperature.Set (n8, v8);
  m_temperature.Set (n9, v9);
}

void
ReliabilityHelper::SetApplication(std::string n0, const DoubleValue &v0)
{
  m_appname = n0;
  m_datasize = v0.Get();
}
void
ReliabilityHelper::Install (Ptr<Node> node)
{
  Ptr<Object> object = node;

  Ptr<PowerModel> powermodel = object->GetObject<PowerModel> ();
  if (powermodel == 0)
    {
      powermodel = m_power.Create ()->GetObject<PowerModel> ();
      if (powermodel == 0)
        {
          NS_FATAL_ERROR ("The requested power model is not a valid power model: \""<< 
                          m_power.GetTypeId ().GetName ()<<"\"");
        }
      NS_LOG_DEBUG ("node="<<object<<", mob="<<powermodel);
      object->AggregateObject (powermodel);
    }

  Ptr<PerformanceModel> performancemodel = object->GetObject<PerformanceModel> ();
  if (performancemodel == 0)
    {
      performancemodel = m_performance.Create ()->GetObject<PerformanceModel> ();
      if (performancemodel == 0)
        {
          NS_FATAL_ERROR ("The requested power model is not a valid performance model: \""<< 
                          m_performance.GetTypeId ().GetName ()<<"\"");
        }
      NS_LOG_DEBUG ("node="<<object<<", mob="<<performancemodel);
      object->AggregateObject (performancemodel);
    }

  Ptr<TemperatureModel> temperaturemodel = object->GetObject<TemperatureModel> ();
  if (temperaturemodel == 0)
    {
      temperaturemodel = m_temperature.Create ()->GetObject<TemperatureModel> ();
      if (temperaturemodel == 0)
        {
          NS_FATAL_ERROR ("The requested power model is not a valid temperature model: \""<< 
                          m_temperature.GetTypeId ().GetName ()<<"\"");
        }
      NS_LOG_DEBUG ("node="<<object<<", mob="<<temperaturemodel);
      object->AggregateObject (temperaturemodel);
    }

    powermodel->RegisterPerformanceModel(performancemodel);
    powermodel->RegisterTemperatureModel (temperaturemodel);
    powermodel->SetApplication(m_appname,m_datasize);
    
}

} // namespace ns3
