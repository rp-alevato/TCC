/***************************************************************************//**
 * @file
 * @brief board_features.h
 *******************************************************************************
 * # License
 * <b>Copyright 2018 Silicon Laboratories Inc. www.silabs.com</b>
 *******************************************************************************
 *
 * The licensor of this software is Silicon Laboratories Inc. Your use of this
 * software is governed by the terms of Silicon Labs Master Software License
 * Agreement (MSLA) available at
 * www.silabs.com/about-us/legal/master-software-license-agreement. This
 * software is distributed to you in Source Code format and is governed by the
 * sections of the MSLA applicable to Source Code.
 *
 ******************************************************************************/

#ifndef BOARD_FEATURES_H
#define BOARD_FEATURES_H

#include "ble-configuration.h"

/* Indicate if LCD is supported */
#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4171A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180B)\
  || (EMBER_AF_BOARD_TYPE == BRD4181A)\
  || (EMBER_AF_BOARD_TYPE == BRD4181B)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_LCD_SUPPORT
#endif

/* Indicate if the same pins are used for LEDs and Buttons on the WSTK */
#if (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4309A)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_LED_BUTTON_ON_SAME_PIN
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4160A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4166A)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4171A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183B)\
  || (EMBER_AF_BOARD_TYPE == BRD4184A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_SPI_FLASH
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4101B)
#define FEATURE_IOEXPANDER
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4171A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4180A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180B)\
  || (EMBER_AF_BOARD_TYPE == BRD4181A)\
  || (EMBER_AF_BOARD_TYPE == BRD4181B)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4309A)
#define FEATURE_PA_INPUT_FROM_VBAT
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183B)\
  || (EMBER_AF_BOARD_TYPE == BRD4184A)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_EXP_HEADER_USART1
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)
#define FEATURE_EXP_HEADER_USART3
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4160A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4166A)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4171A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4179B)\
  || (EMBER_AF_BOARD_TYPE == BRD4180A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180B)\
  || (EMBER_AF_BOARD_TYPE == BRD4181A)\
  || (EMBER_AF_BOARD_TYPE == BRD4181B)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183B)\
  || (EMBER_AF_BOARD_TYPE == BRD4184A)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4309A)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_PTI_SUPPORT
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4179B)\
  || (EMBER_AF_BOARD_TYPE == BRD4180A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180B)\
  || (EMBER_AF_BOARD_TYPE == BRD4181A)\
  || (EMBER_AF_BOARD_TYPE == BRD4181B)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183B)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4309A)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_HW_FLOW_CONTROL
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4160A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4166A)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4171A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4184A)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)
#define FEATURE_I2C_SENSOR
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4160A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4166A)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4179B)\
  || (EMBER_AF_BOARD_TYPE == BRD4180A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180B)\
  || (EMBER_AF_BOARD_TYPE == BRD4181A)\
  || (EMBER_AF_BOARD_TYPE == BRD4181B)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183B)\
  || (EMBER_AF_BOARD_TYPE == BRD4184A)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)
#define FEATURE_LFXO
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4179B)
#define FEATURE_EFP
#endif

#if (EMBER_AF_BOARD_TYPE == BRD4100A)\
  || (EMBER_AF_BOARD_TYPE == BRD4101B)\
  || (EMBER_AF_BOARD_TYPE == BRD4103A)\
  || (EMBER_AF_BOARD_TYPE == BRD4104A)\
  || (EMBER_AF_BOARD_TYPE == BRD4105A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150A)\
  || (EMBER_AF_BOARD_TYPE == BRD4150B)\
  || (EMBER_AF_BOARD_TYPE == BRD4150C)\
  || (EMBER_AF_BOARD_TYPE == BRD4151A)\
  || (EMBER_AF_BOARD_TYPE == BRD4153A)\
  || (EMBER_AF_BOARD_TYPE == BRD4158A)\
  || (EMBER_AF_BOARD_TYPE == BRD4159A)\
  || (EMBER_AF_BOARD_TYPE == BRD4160A)\
  || (EMBER_AF_BOARD_TYPE == BRD4161A)\
  || (EMBER_AF_BOARD_TYPE == BRD4162A)\
  || (EMBER_AF_BOARD_TYPE == BRD4163A)\
  || (EMBER_AF_BOARD_TYPE == BRD4164A)\
  || (EMBER_AF_BOARD_TYPE == BRD4165B)\
  || (EMBER_AF_BOARD_TYPE == BRD4166A)\
  || (EMBER_AF_BOARD_TYPE == BRD4167A)\
  || (EMBER_AF_BOARD_TYPE == BRD4168A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169A)\
  || (EMBER_AF_BOARD_TYPE == BRD4169B)\
  || (EMBER_AF_BOARD_TYPE == BRD4170A)\
  || (EMBER_AF_BOARD_TYPE == BRD4171A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172A)\
  || (EMBER_AF_BOARD_TYPE == BRD4172B)\
  || (EMBER_AF_BOARD_TYPE == BRD4173A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174A)\
  || (EMBER_AF_BOARD_TYPE == BRD4174B)\
  || (EMBER_AF_BOARD_TYPE == BRD4175A)\
  || (EMBER_AF_BOARD_TYPE == BRD4179B)\
  || (EMBER_AF_BOARD_TYPE == BRD4180A)\
  || (EMBER_AF_BOARD_TYPE == BRD4180B)\
  || (EMBER_AF_BOARD_TYPE == BRD4181A)\
  || (EMBER_AF_BOARD_TYPE == BRD4181B)\
  || (EMBER_AF_BOARD_TYPE == BRD4182A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183A)\
  || (EMBER_AF_BOARD_TYPE == BRD4183B)\
  || (EMBER_AF_BOARD_TYPE == BRD4184A)\
  || (EMBER_AF_BOARD_TYPE == BRD4300A)\
  || (EMBER_AF_BOARD_TYPE == BRD4301A)\
  || (EMBER_AF_BOARD_TYPE == BRD4302A)\
  || (EMBER_AF_BOARD_TYPE == BRD4303A)\
  || (EMBER_AF_BOARD_TYPE == BRD4304A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305A)\
  || (EMBER_AF_BOARD_TYPE == BRD4305C)\
  || (EMBER_AF_BOARD_TYPE == BRD4305D)\
  || (EMBER_AF_BOARD_TYPE == BRD4305E)\
  || (EMBER_AF_BOARD_TYPE == BRD4306A)\
  || (EMBER_AF_BOARD_TYPE == BRD4306B)\
  || (EMBER_AF_BOARD_TYPE == BRD4306C)\
  || (EMBER_AF_BOARD_TYPE == BRD4306D)\
  || (EMBER_AF_BOARD_TYPE == BRD4308A)\
  || (EMBER_AF_BOARD_TYPE == BRD4308B)\
  || (EMBER_AF_BOARD_TYPE == BRD4309A)\
  || (EMBER_AF_BOARD_TYPE == BRD4310A)\
  || (EMBER_AF_BOARD_TYPE == BRD4311A)\
  || (EMBER_AF_BOARD_TYPE == RD_0057_0101)
#define FEATURE_BOARD_DETECTED
#endif

#if (EMBER_AF_BOARD_TYPE == CUSTOM_BOARD)
// Uncomment the corresponding line in case of using Silicon Labs board feature in your design.
// For using the selected feature you may need additional drivers. Check an appropriate SDK example for reference.

// #define FEATURE_LCD_SUPPORT
// #define FEATURE_LED_BUTTON_ON_SAME_PIN
// #define FEATURE_SPI_FLASH
// #define FEATURE_IOEXPANDER
// #define FEATURE_PA_INPUT_FROM_VBAT
// #define FEATURE_EXP_HEADER_USART1
// #define FEATURE_EXP_HEADER_USART3
// #define FEATURE_PTI_SUPPORT
// #define FEATURE_HW_FLOW_CONTROL
// #define FEATURE_I2C_SENSOR
// #define FEATURE_LFXO
// #define FEATURE_EFP
// #define FEATURE_BOARD_DETECTED
#endif

#endif /* BOARD_FEATURES_H */