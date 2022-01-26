/***************************************************************************//**
 * @file
 * @brief Silabs Network Co-Processor (NCP) Empty Target Sample Application.
 * This application allows the user to use as an NCP target device if connected
 * to an appropriate host application.
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

#include "init_mcu.h"
#include "init_board.h"
#include "init_app.h"
#include "ble-configuration.h"
#include "board_features.h"

/* BG stack headers */
#include "bg_types.h"
#include "ncp_gecko.h"
#include "gatt_db.h"
#include "ncp_usart.h"
#include "em_core.h"

/* libraries containing default gecko configuration values */
#include "em_emu.h"
#include "em_cmu.h"

/* Device initialization header */
#include "hal-config.h"

#ifdef FEATURE_BOARD_DETECTED
#if defined(HAL_CONFIG)
#include "bsphalconfig.h"
#else
#include "bspconfig.h"
#endif // HAL_CONFIG
#endif // FEATURE_BOARD_DETECTED

/***********************************************************************************************//**
 * @addtogroup Application
 * @{
 **************************************************************************************************/

/***********************************************************************************************//**
 * @addtogroup app
 * @{
 **************************************************************************************************/

#ifndef MAX_CONNECTIONS
#define MAX_CONNECTIONS 4
#endif
#ifndef MAX_ADVERTISERS
#define MAX_ADVERTISERS 4
#endif
#ifndef MAX_PERIODIC_SYNC
#define MAX_PERIODIC_SYNC 4
#endif

#if defined(NCP_DEEP_SLEEP_ENABLED)
#define NCP_TX_WATCH_COUNTER_MAX    1000
#else
#define NCP_TX_WATCH_COUNTER_MAX    0xffff
#endif

uint8_t bluetooth_stack_heap[DEFAULT_BLUETOOTH_HEAP(MAX_CONNECTIONS + MAX_PERIODIC_SYNC)];
static int ncp_tx_queue_watch_counter = 0;

// Gecko configuration parameters (see gecko_configuration.h)
static const gecko_configuration_t config = {
  .config_flags = 0,
#if (defined(FEATURE_LFXO) || defined(PLFRCO_PRESENT)) && defined(NCP_DEEP_SLEEP_ENABLED)
  .sleep.flags = SLEEP_FLAGS_DEEP_SLEEP_ENABLE,
#else
  .sleep.flags = 0,
#endif
  .bluetooth.max_connections = MAX_CONNECTIONS,
  .bluetooth.max_advertisers = MAX_ADVERTISERS,
  .bluetooth.max_periodic_sync = MAX_PERIODIC_SYNC,
  .bluetooth.heap = bluetooth_stack_heap,
  .bluetooth.heap_size = sizeof(bluetooth_stack_heap),
#if defined(FEATURE_LFXO)
  .bluetooth.sleep_clock_accuracy = 100, // ppm
#elif defined(PLFRCO_PRESENT)
  .bluetooth.sleep_clock_accuracy = 500, // ppm
#endif
  .gattdb = &bg_gattdb_data,
  .pa.config_enable = 1, // Set this to be a valid PA config
#if defined(FEATURE_PA_INPUT_FROM_VBAT)
  .pa.input = GECKO_RADIO_PA_INPUT_VBAT, // Configure PA input to VBAT
#else
  .pa.input = GECKO_RADIO_PA_INPUT_DCDC,
#endif // defined(FEATURE_PA_INPUT_FROM_VBAT)
  .rf.flags = GECKO_RF_CONFIG_ANTENNA,                 /* Enable antenna configuration. */
  .rf.antenna = GECKO_RF_ANTENNA,                      /* Select antenna path! */
};

/**
 * Handle events meant to be handled locally
 */
static uint32_t local_handle_event(struct gecko_cmd_packet *evt)
{
  bool evt_handled = false;
  switch (BGLIB_MSG_ID(evt->header)) {
    default:
      break;
  }
  return evt_handled;
}

/*
 * Check if the TX queue has enough space for a new event message.
 * A counter is used here for catching a rare situation where the NCP gets stuck somehow.
 */
static bool check_ncp_tx_queue_space()
{
  if (ncp_get_transmit_space_for_events() < NCP_CMD_SIZE) {
    ++ncp_tx_queue_watch_counter;
    if (ncp_tx_queue_watch_counter > NCP_TX_WATCH_COUNTER_MAX) {
      // Fatal error, or the NCP transport has been shut down for too long.
      // Reset the queue for keeping the device responsive in case the stack is still running in background.
      tx_queue_reset();
      ncp_tx_queue_watch_counter = 0;
      return true;
    }
    return false;
  } else {
    ncp_tx_queue_watch_counter = 0;
    return true;
  }
}

/**
 * Extract the next BGAPI event only when there is enough space for buffering it in TX queue.
 */
static struct gecko_cmd_packet *get_bgapi_event()
{
  if (check_ncp_tx_queue_space()) {
    return gecko_peek_event();
  }
  return NULL;
}

int main(void)
{
  // Initialize device
  initMcu();
  // Initialize board
  initBoard();
  // Initialize application
  initApp();

  // Initialize stack
  gecko_init(&config);

  // Initialize periodic advertisement scanner in Bluetooth stack
  ncp_gecko_bgapi_class_sync_init();

  // Initialize CTE receiver class in Bluetooth stack
  ncp_gecko_bgapi_class_cte_receiver_init();

  // NCP USART init
  ncp_usart_init();
  initVcomEnable();

  while (1) {
    struct gecko_cmd_packet *evt;
    ncp_handle_command();
    for (evt = get_bgapi_event(); evt != NULL; evt = get_bgapi_event()) {
      if (!ncp_handle_event(evt) && !local_handle_event(evt)) {
        // send out the event if not handled either by NCP or locally
        // filter out scan responses to reduce uart traffic
        if (BGLIB_MSG_ID(evt->header) != gecko_evt_le_gap_scan_response_id) {
          ncp_transmit_enqueue(evt);
        }
      }
      // if a command is received, break and handle the command
      if (ncp_command_received()) {
        break;
      }
    }
    ncp_transmit();

    // If an NCP command received, do not sleep
    if (!ncp_command_received()) {
      CORE_CRITICAL_IRQ_DISABLE();
      gecko_sleep_for_ms(gecko_can_sleep_ms());
      CORE_CRITICAL_IRQ_ENABLE();
    }
  }
}

/** @} (end addtogroup app) */
/** @} (end addtogroup Application) */
