/***************************************************************************//**
 * @brief Bluetooth BGAPI for host applications in NCP mode
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

#ifndef HOST_GECKO_H
#define HOST_GECKO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include "bg_types.h"
#include "gecko_configuration.h"
#include "bg_errorcodes.h"
typedef uint8_t   uint8;
typedef uint16_t  uint16;
typedef uint32_t  uint32;
typedef uint64_t  uint64;
typedef int8_t    int8;
typedef int16_t   int16;
typedef int32_t   int32;
typedef int64_t   int64;

typedef struct {
  uint16_t len;
  uint8_t  data[];
}uint16array;



/* Compatibility */
#ifndef PACKSTRUCT
/*Default packed configuration*/
#ifdef __GNUC__
#ifdef _WIN32
#define PACKSTRUCT( decl ) decl __attribute__((__packed__,gcc_struct))
#else
#define PACKSTRUCT( decl ) decl __attribute__((__packed__))
#endif
#define ALIGNED __attribute__((aligned(0x4)))
#elif __IAR_SYSTEMS_ICC__

#define PACKSTRUCT( decl ) __packed decl

#define ALIGNED
#elif _MSC_VER  /*msvc*/

#define PACKSTRUCT( decl ) __pragma( pack(push, 1) ) decl __pragma( pack(pop) )
#define ALIGNED
#else 
#define PACKSTRUCT(a) a PACKED 
#endif
#endif


#define BGLIB_DEPRECATED_API __attribute__((deprecated))
#define BGLIB_MSG_ID(HDR) ((HDR)&0xffff00f8)
#define BGLIB_MSG_HEADER_LEN (4)
#define BGLIB_MSG_LEN(HDR) ((((HDR)&0x7)<<8)|(((HDR)&0xff00)>>8))

/**
 * The maximum BGAPI command payload size.
 */
#define BGLIB_MSG_MAX_PAYLOAD 256


#define BGLIB_BIT_ENCRYPTED (1 << 6) // Bit indicating whether the packet is encrypted
#define BGLIB_MSG_ENCRYPTED(HDR) ((HDR) & BGLIB_BIT_ENCRYPTED)
/**
 * Blocks until new event arrives which requires processing by user application.
 * 
 * @return pointer to received event
 */
struct gecko_cmd_packet* gecko_wait_event(void);

/**
 * Same as gecko_wait_event but does not block if no events waiting, instead returns NULL
 *
 * @return pointer to received event or NULL if no event waiting
 */
struct gecko_cmd_packet* gecko_peek_event(void);

/**
 * Events are in queue waiting for processing
 * Call gecko_wait_event or gecko_peek_event to process pending events
 *
 * @return nonzero if processing required
 */
int gecko_event_pending(void);

/**
 * Initialize stack.
 * @param config The pointer to configuration parameters, cannot be NULL.
 * @return bg_err_success if the initialization was successful; Other error code
 *         indicates a failure on initializing persistent storage.
 */
errorcode_t gecko_stack_init(const gecko_configuration_t *config);







enum system_linklayer_config_key
{
    system_linklayer_config_key_halt                             = 0x1,
    system_linklayer_config_key_priority_range                   = 0x2,
    system_linklayer_config_key_scan_channels                    = 0x3,
    system_linklayer_config_key_set_flags                        = 0x4,
    system_linklayer_config_key_clr_flags                        = 0x5,
    system_linklayer_config_key_set_afh_interval                 = 0x7,
    system_linklayer_config_key_set_priority_table               = 0x9
};

enum le_gap_address_type
{
    le_gap_address_type_public                                   = 0x0,
    le_gap_address_type_random                                   = 0x1
};

enum le_gap_phy_type
{
    le_gap_phy_1m                                                = 0x1,
    le_gap_phy_2m                                                = 0x2,
    le_gap_phy_coded                                             = 0x4
};

enum le_gap_connectable_mode
{
    le_gap_non_connectable                                       = 0x0,
    le_gap_directed_connectable                                  = 0x1,
    le_gap_undirected_connectable                                = 0x2,
    le_gap_connectable_scannable                                 = 0x2,
    le_gap_scannable_non_connectable                             = 0x3,
    le_gap_connectable_non_scannable                             = 0x4
};

enum le_gap_discoverable_mode
{
    le_gap_non_discoverable                                      = 0x0,
    le_gap_limited_discoverable                                  = 0x1,
    le_gap_general_discoverable                                  = 0x2,
    le_gap_broadcast                                             = 0x3,
    le_gap_user_data                                             = 0x4
};

enum le_gap_discover_mode
{
    le_gap_discover_limited                                      = 0x0,
    le_gap_discover_generic                                      = 0x1,
    le_gap_discover_observation                                  = 0x2
};

enum le_gap_adv_address_type
{
    le_gap_identity_address                                      = 0x0,
    le_gap_non_resolvable                                        = 0x1
};

enum sync_advertiser_clock_accuracy
{
    sync_clock_accuracy_500                                      = 0x1f4,
    sync_clock_accuracy_250                                      = 0xfa,
    sync_clock_accuracy_150                                      = 0x96,
    sync_clock_accuracy_100                                      = 0x64,
    sync_clock_accuracy_75                                       = 0x4b,
    sync_clock_accuracy_50                                       = 0x32,
    sync_clock_accuracy_30                                       = 0x1e,
    sync_clock_accuracy_20                                       = 0x14
};

enum le_connection_security
{
    le_connection_mode1_level1                                   = 0x0,
    le_connection_mode1_level2                                   = 0x1,
    le_connection_mode1_level3                                   = 0x2,
    le_connection_mode1_level4                                   = 0x3
};

enum gatt_att_opcode
{
    gatt_read_by_type_request                                    = 0x8,
    gatt_read_by_type_response                                   = 0x9,
    gatt_read_request                                            = 0xa,
    gatt_read_response                                           = 0xb,
    gatt_read_blob_request                                       = 0xc,
    gatt_read_blob_response                                      = 0xd,
    gatt_read_multiple_request                                   = 0xe,
    gatt_read_multiple_response                                  = 0xf,
    gatt_write_request                                           = 0x12,
    gatt_write_response                                          = 0x13,
    gatt_write_command                                           = 0x52,
    gatt_prepare_write_request                                   = 0x16,
    gatt_prepare_write_response                                  = 0x17,
    gatt_execute_write_request                                   = 0x18,
    gatt_execute_write_response                                  = 0x19,
    gatt_handle_value_notification                               = 0x1b,
    gatt_handle_value_indication                                 = 0x1d
};

enum gatt_client_config_flag
{
    gatt_disable                                                 = 0x0,
    gatt_notification                                            = 0x1,
    gatt_indication                                              = 0x2
};

enum gatt_execute_write_flag
{
    gatt_cancel                                                  = 0x0,
    gatt_commit                                                  = 0x1
};

enum gatt_server_characteristic_status_flag
{
    gatt_server_client_config                                    = 0x1,
    gatt_server_confirmation                                     = 0x2
};


enum test_packet_type
{
    test_pkt_prbs9                                               = 0x0,
    test_pkt_11110000                                            = 0x1,
    test_pkt_10101010                                            = 0x2,
    test_pkt_11111111                                            = 0x4,
    test_pkt_00000000                                            = 0x5,
    test_pkt_00001111                                            = 0x6,
    test_pkt_01010101                                            = 0x7,
    test_pkt_pn9                                                 = 0xfd,
    test_pkt_carrier                                             = 0xfe
};

enum test_phy
{
    test_phy_1m                                                  = 0x1,
    test_phy_2m                                                  = 0x2,
    test_phy_125k                                                = 0x3,
    test_phy_500k                                                = 0x4
};

enum sm_bonding_key
{
    sm_bonding_key_ltk                                           = 0x1,
    sm_bonding_key_addr_public                                   = 0x2,
    sm_bonding_key_addr_static                                   = 0x4,
    sm_bonding_key_irk                                           = 0x8,
    sm_bonding_key_edivrand                                      = 0x10,
    sm_bonding_key_csrk                                          = 0x20,
    sm_bonding_key_masterid                                      = 0x40
};

enum sm_io_capability
{
    sm_io_capability_displayonly                                 = 0x0,
    sm_io_capability_displayyesno                                = 0x1,
    sm_io_capability_keyboardonly                                = 0x2,
    sm_io_capability_noinputnooutput                             = 0x3,
    sm_io_capability_keyboarddisplay                             = 0x4
};

enum homekit_category
{
    homekit_not_allowed                                          = 0x0,
    homekit_other                                                = 0x1,
    homekit_bridge                                               = 0x2,
    homekit_fan                                                  = 0x3,
    homekit_garage                                               = 0x4,
    homekit_lightbulb                                            = 0x5,
    homekit_doorlock                                             = 0x6,
    homekit_outlet                                               = 0x7,
    homekit_switch_accessory                                     = 0x8,
    homekit_thermostat                                           = 0x9,
    homekit_sensor                                               = 0xa,
    homekit_security_system                                      = 0xb,
    homekit_door                                                 = 0xc,
    homekit_window                                               = 0xd,
    homekit_window_covering                                      = 0xe,
    homekit_programmable_switch                                  = 0xf,
    homekit_ip_camera                                            = 0x11,
    homekit_video_door_bell                                      = 0x12,
    homekit_air_purifier                                         = 0x13,
    homekit_heater                                               = 0x14,
    homekit_air_conditioner                                      = 0x15,
    homekit_humidifier                                           = 0x16,
    homekit_dehumidifier                                         = 0x17,
    homekit_sprinkler                                            = 0x1c,
    homekit_faucet                                               = 0x1d,
    homekit_shower_system                                        = 0x1e,
    homekit_remote                                               = 0x20,
    homekit_wifi_router                                          = 0x21
};

enum homekit_status_code
{
    homekit_success                                              = 0x0,
    homekit_invalid_request                                      = 0x6
};

enum coex_option
{
    coex_option_enable                                           = 0x100,
    coex_option_tx_abort                                         = 0x400,
    coex_option_high_priority                                    = 0x800
};

enum l2cap_coc_connection_result
{
    l2cap_connection_successful                                  = 0x0,
    l2cap_le_psm_not_supported                                   = 0x2,
    l2cap_no_resources_available                                 = 0x4,
    l2cap_insufficient_authentication                            = 0x5,
    l2cap_insufficient_authorization                             = 0x6,
    l2cap_insufficient_encryption_key_size                       = 0x7,
    l2cap_insufficient_encryption                                = 0x8,
    l2cap_invalid_source_cid                                     = 0x9,
    l2cap_source_cid_already_allocated                           = 0xa,
    l2cap_unacceptable_parameters                                = 0xb
};

enum l2cap_command_reject_reason
{
    l2cap_command_not_understood                                 = 0x0,
    l2cap_signaling_mtu_exceeded                                 = 0x1,
    l2cap_invalid_cid_request                                    = 0x2
};

enum l2cap_command_code
{
    l2cap_disconnection_request                                  = 0x6,
    l2cap_connection_request                                     = 0x14,
    l2cap_flow_control_credit                                    = 0x16
};

enum gecko_parameter_types
{
    gecko_msg_parameter_uint8=2,
    gecko_msg_parameter_int8=3,
    gecko_msg_parameter_uint16=4,
    gecko_msg_parameter_int16=5,
    gecko_msg_parameter_uint32=6,
    gecko_msg_parameter_int32=7,
    gecko_msg_parameter_uint8array=8,
    gecko_msg_parameter_string=9,
    gecko_msg_parameter_hwaddr=10,
    gecko_msg_parameter_uint16array=11
};

enum gecko_msg_types
{
    gecko_msg_type_cmd=0x00,
    gecko_msg_type_rsp=0x00,
    gecko_msg_type_evt=0x80
};
enum gecko_dev_types
{
    gecko_dev_type_gecko   =0x20
};



#define FLASH_PS_KEY_LOCAL_BD_ADDR                                   0x2c
#define FLASH_PS_KEY_TX_POWER                                        0x31
#define FLASH_PS_KEY_CTUNE                                           0x32
#define FLASH_PS_KEY_APPLICATION_GSN                                 0x33
#define FLASH_PS_KEY_OTA_FLAGS                                       0x35
#define FLASH_PS_KEY_OTA_DEVICE_NAME                                 0x36
#define FLASH_PS_KEY_DEVICE_IRK                                      0x37
#define FLASH_PS_KEY_BONDING_PRIORITY_LIST                           0x38
#define FLASH_PS_KEY_OTA_ADVERTISEMENT_PACKET                        0x39
#define FLASH_PS_KEY_OTA_SCAN_RESPONSE_PACKET                        0x3a
#define FLASH_PS_KEY_APPLICATION_AI                                  0x3b
#define FLASH_PS_KEY_IDENTITY_ADDR_TYPE                              0x3c
#define FLASH_PS_KEY_GATT_DB_HASH                                    0x3d
#define FLASH_PS_KEY_OTA_RF_PATH                                     0x3e
#define FLASH_PS_KEY_BONDING_DB_CONFIG                               0x3fff


#define gecko_cmd_dfu_reset_id                                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00000000)
#define gecko_cmd_dfu_flash_set_address_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01000000)
#define gecko_cmd_dfu_flash_upload_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02000000)
#define gecko_cmd_dfu_flash_upload_finish_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03000000)
#define gecko_cmd_system_hello_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00010000)
#define gecko_cmd_system_reset_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01010000)
#define gecko_cmd_system_get_bt_address_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03010000)
#define gecko_cmd_system_set_bt_address_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04010000)
#define gecko_cmd_system_set_tx_power_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0a010000)
#define gecko_cmd_system_get_random_data_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0b010000)
#define gecko_cmd_system_halt_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0c010000)
#define gecko_cmd_system_set_device_name_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0d010000)
#define gecko_cmd_system_linklayer_configure_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0e010000)
#define gecko_cmd_system_get_counters_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0f010000)
#define gecko_cmd_system_data_buffer_write_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x12010000)
#define gecko_cmd_system_set_identity_address_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x13010000)
#define gecko_cmd_system_data_buffer_clear_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x14010000)
#define gecko_cmd_le_gap_open_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00030000)
#define gecko_cmd_le_gap_set_mode_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01030000)
#define gecko_cmd_le_gap_discover_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02030000)
#define gecko_cmd_le_gap_end_procedure_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03030000)
#define gecko_cmd_le_gap_set_adv_parameters_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04030000)
#define gecko_cmd_le_gap_set_conn_parameters_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05030000)
#define gecko_cmd_le_gap_set_scan_parameters_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x06030000)
#define gecko_cmd_le_gap_set_adv_data_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x07030000)
#define gecko_cmd_le_gap_set_adv_timeout_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x08030000)
#define gecko_cmd_le_gap_set_conn_phy_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x09030000)
#define gecko_cmd_le_gap_bt5_set_mode_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0a030000)
#define gecko_cmd_le_gap_bt5_set_adv_parameters_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0b030000)
#define gecko_cmd_le_gap_bt5_set_adv_data_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0c030000)
#define gecko_cmd_le_gap_set_privacy_mode_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0d030000)
#define gecko_cmd_le_gap_set_advertise_timing_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0e030000)
#define gecko_cmd_le_gap_set_advertise_channel_map_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0f030000)
#define gecko_cmd_le_gap_set_advertise_report_scan_request_id         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x10030000)
#define gecko_cmd_le_gap_set_advertise_phy_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x11030000)
#define gecko_cmd_le_gap_set_advertise_configuration_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x12030000)
#define gecko_cmd_le_gap_clear_advertise_configuration_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x13030000)
#define gecko_cmd_le_gap_start_advertising_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x14030000)
#define gecko_cmd_le_gap_stop_advertising_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x15030000)
#define gecko_cmd_le_gap_set_discovery_timing_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x16030000)
#define gecko_cmd_le_gap_set_discovery_type_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x17030000)
#define gecko_cmd_le_gap_start_discovery_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x18030000)
#define gecko_cmd_le_gap_set_data_channel_classification_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x19030000)
#define gecko_cmd_le_gap_connect_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x1a030000)
#define gecko_cmd_le_gap_set_advertise_tx_power_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x1b030000)
#define gecko_cmd_le_gap_set_discovery_extended_scan_response_id      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x1c030000)
#define gecko_cmd_le_gap_start_periodic_advertising_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x1d030000)
#define gecko_cmd_le_gap_stop_periodic_advertising_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x1f030000)
#define gecko_cmd_le_gap_set_long_advertising_data_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x20030000)
#define gecko_cmd_le_gap_enable_whitelisting_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x21030000)
#define gecko_cmd_le_gap_set_conn_timing_parameters_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x22030000)
#define gecko_cmd_le_gap_set_advertise_random_address_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x25030000)
#define gecko_cmd_le_gap_clear_advertise_random_address_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x26030000)
#define gecko_cmd_sync_open_id                                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00420000)
#define gecko_cmd_sync_close_id                                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01420000)
#define gecko_cmd_le_connection_set_parameters_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00080000)
#define gecko_cmd_le_connection_get_rssi_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01080000)
#define gecko_cmd_le_connection_disable_slave_latency_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02080000)
#define gecko_cmd_le_connection_set_phy_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03080000)
#define gecko_cmd_le_connection_close_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04080000)
#define gecko_cmd_le_connection_set_timing_parameters_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05080000)
#define gecko_cmd_le_connection_read_channel_map_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x06080000)
#define gecko_cmd_le_connection_set_preferred_phy_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x07080000)
#define gecko_cmd_gatt_set_max_mtu_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00090000)
#define gecko_cmd_gatt_discover_primary_services_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01090000)
#define gecko_cmd_gatt_discover_primary_services_by_uuid_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02090000)
#define gecko_cmd_gatt_discover_characteristics_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03090000)
#define gecko_cmd_gatt_discover_characteristics_by_uuid_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04090000)
#define gecko_cmd_gatt_set_characteristic_notification_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05090000)
#define gecko_cmd_gatt_discover_descriptors_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x06090000)
#define gecko_cmd_gatt_read_characteristic_value_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x07090000)
#define gecko_cmd_gatt_read_characteristic_value_by_uuid_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x08090000)
#define gecko_cmd_gatt_write_characteristic_value_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x09090000)
#define gecko_cmd_gatt_write_characteristic_value_without_response_id  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0a090000)
#define gecko_cmd_gatt_prepare_characteristic_value_write_id          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0b090000)
#define gecko_cmd_gatt_execute_characteristic_value_write_id          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0c090000)
#define gecko_cmd_gatt_send_characteristic_confirmation_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0d090000)
#define gecko_cmd_gatt_read_descriptor_value_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0e090000)
#define gecko_cmd_gatt_write_descriptor_value_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0f090000)
#define gecko_cmd_gatt_find_included_services_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x10090000)
#define gecko_cmd_gatt_read_multiple_characteristic_values_id         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x11090000)
#define gecko_cmd_gatt_read_characteristic_value_from_offset_id       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x12090000)
#define gecko_cmd_gatt_prepare_characteristic_value_reliable_write_id  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x13090000)
#define gecko_cmd_gatt_server_read_attribute_value_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x000a0000)
#define gecko_cmd_gatt_server_read_attribute_type_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x010a0000)
#define gecko_cmd_gatt_server_write_attribute_value_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x020a0000)
#define gecko_cmd_gatt_server_send_user_read_response_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x030a0000)
#define gecko_cmd_gatt_server_send_user_write_response_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x040a0000)
#define gecko_cmd_gatt_server_send_characteristic_notification_id     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x050a0000)
#define gecko_cmd_gatt_server_find_attribute_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x060a0000)
#define gecko_cmd_gatt_server_set_capabilities_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x080a0000)
#define gecko_cmd_gatt_server_set_max_mtu_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0a0a0000)
#define gecko_cmd_gatt_server_get_mtu_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0b0a0000)
#define gecko_cmd_gatt_server_enable_capabilities_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0c0a0000)
#define gecko_cmd_gatt_server_disable_capabilities_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0d0a0000)
#define gecko_cmd_gatt_server_get_enabled_capabilities_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0e0a0000)
#define gecko_cmd_hardware_set_soft_timer_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x000c0000)
#define gecko_cmd_hardware_get_time_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0b0c0000)
#define gecko_cmd_hardware_set_lazy_soft_timer_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0c0c0000)
#define gecko_cmd_flash_ps_erase_all_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x010d0000)
#define gecko_cmd_flash_ps_save_id                                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x020d0000)
#define gecko_cmd_flash_ps_load_id                                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x030d0000)
#define gecko_cmd_flash_ps_erase_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x040d0000)
#define gecko_cmd_test_dtm_tx_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x000e0000)
#define gecko_cmd_test_dtm_rx_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x010e0000)
#define gecko_cmd_test_dtm_end_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x020e0000)
#define gecko_cmd_sm_set_bondable_mode_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x000f0000)
#define gecko_cmd_sm_configure_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x010f0000)
#define gecko_cmd_sm_store_bonding_configuration_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x020f0000)
#define gecko_cmd_sm_increase_security_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x040f0000)
#define gecko_cmd_sm_delete_bonding_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x060f0000)
#define gecko_cmd_sm_delete_bondings_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x070f0000)
#define gecko_cmd_sm_enter_passkey_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x080f0000)
#define gecko_cmd_sm_passkey_confirm_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x090f0000)
#define gecko_cmd_sm_set_oob_data_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0a0f0000)
#define gecko_cmd_sm_list_all_bondings_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0b0f0000)
#define gecko_cmd_sm_bonding_confirm_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0e0f0000)
#define gecko_cmd_sm_set_debug_mode_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0f0f0000)
#define gecko_cmd_sm_set_passkey_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x100f0000)
#define gecko_cmd_sm_use_sc_oob_id                                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x110f0000)
#define gecko_cmd_sm_set_sc_remote_oob_data_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x120f0000)
#define gecko_cmd_sm_add_to_whitelist_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x130f0000)
#define gecko_cmd_sm_set_minimum_key_size_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x140f0000)
#define gecko_cmd_homekit_configure_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00130000)
#define gecko_cmd_homekit_advertise_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01130000)
#define gecko_cmd_homekit_delete_pairings_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02130000)
#define gecko_cmd_homekit_check_authcp_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03130000)
#define gecko_cmd_homekit_get_pairing_id_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04130000)
#define gecko_cmd_homekit_send_write_response_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05130000)
#define gecko_cmd_homekit_send_read_response_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x06130000)
#define gecko_cmd_homekit_gsn_action_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x07130000)
#define gecko_cmd_homekit_event_notification_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x08130000)
#define gecko_cmd_homekit_broadcast_action_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x09130000)
#define gecko_cmd_homekit_configure_product_data_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x0a130000)
#define gecko_cmd_coex_set_options_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00200000)
#define gecko_cmd_coex_get_counters_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01200000)
#define gecko_cmd_coex_set_parameters_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02200000)
#define gecko_cmd_coex_set_directional_priority_pulse_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03200000)
#define gecko_cmd_l2cap_coc_send_connection_request_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01430000)
#define gecko_cmd_l2cap_coc_send_connection_response_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02430000)
#define gecko_cmd_l2cap_coc_send_le_flow_control_credit_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03430000)
#define gecko_cmd_l2cap_coc_send_disconnection_request_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04430000)
#define gecko_cmd_l2cap_coc_send_data_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05430000)
#define gecko_cmd_cte_transmitter_enable_connection_cte_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00440000)
#define gecko_cmd_cte_transmitter_disable_connection_cte_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01440000)
#define gecko_cmd_cte_transmitter_enable_connectionless_cte_id        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02440000)
#define gecko_cmd_cte_transmitter_disable_connectionless_cte_id       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03440000)
#define gecko_cmd_cte_transmitter_set_dtm_parameters_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04440000)
#define gecko_cmd_cte_transmitter_clear_dtm_parameters_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05440000)
#define gecko_cmd_cte_transmitter_enable_silabs_cte_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x06440000)
#define gecko_cmd_cte_transmitter_disable_silabs_cte_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x07440000)
#define gecko_cmd_cte_receiver_configure_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00450000)
#define gecko_cmd_cte_receiver_enable_connection_cte_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x01450000)
#define gecko_cmd_cte_receiver_disable_connection_cte_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x02450000)
#define gecko_cmd_cte_receiver_enable_connectionless_cte_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x03450000)
#define gecko_cmd_cte_receiver_disable_connectionless_cte_id          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x04450000)
#define gecko_cmd_cte_receiver_set_dtm_parameters_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x05450000)
#define gecko_cmd_cte_receiver_clear_dtm_parameters_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x06450000)
#define gecko_cmd_cte_receiver_enable_silabs_cte_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x07450000)
#define gecko_cmd_cte_receiver_disable_silabs_cte_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x08450000)
#define gecko_cmd_user_message_to_target_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_cmd|0x00ff0000)

#define gecko_rsp_dfu_reset_id                                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00000000)
#define gecko_rsp_dfu_flash_set_address_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01000000)
#define gecko_rsp_dfu_flash_upload_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02000000)
#define gecko_rsp_dfu_flash_upload_finish_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03000000)
#define gecko_rsp_system_hello_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00010000)
#define gecko_rsp_system_reset_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01010000)
#define gecko_rsp_system_get_bt_address_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03010000)
#define gecko_rsp_system_set_bt_address_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04010000)
#define gecko_rsp_system_set_tx_power_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0a010000)
#define gecko_rsp_system_get_random_data_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0b010000)
#define gecko_rsp_system_halt_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0c010000)
#define gecko_rsp_system_set_device_name_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0d010000)
#define gecko_rsp_system_linklayer_configure_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0e010000)
#define gecko_rsp_system_get_counters_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0f010000)
#define gecko_rsp_system_data_buffer_write_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x12010000)
#define gecko_rsp_system_set_identity_address_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x13010000)
#define gecko_rsp_system_data_buffer_clear_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x14010000)
#define gecko_rsp_le_gap_open_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00030000)
#define gecko_rsp_le_gap_set_mode_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01030000)
#define gecko_rsp_le_gap_discover_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02030000)
#define gecko_rsp_le_gap_end_procedure_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03030000)
#define gecko_rsp_le_gap_set_adv_parameters_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04030000)
#define gecko_rsp_le_gap_set_conn_parameters_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05030000)
#define gecko_rsp_le_gap_set_scan_parameters_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x06030000)
#define gecko_rsp_le_gap_set_adv_data_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x07030000)
#define gecko_rsp_le_gap_set_adv_timeout_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x08030000)
#define gecko_rsp_le_gap_set_conn_phy_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x09030000)
#define gecko_rsp_le_gap_bt5_set_mode_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0a030000)
#define gecko_rsp_le_gap_bt5_set_adv_parameters_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0b030000)
#define gecko_rsp_le_gap_bt5_set_adv_data_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0c030000)
#define gecko_rsp_le_gap_set_privacy_mode_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0d030000)
#define gecko_rsp_le_gap_set_advertise_timing_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0e030000)
#define gecko_rsp_le_gap_set_advertise_channel_map_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0f030000)
#define gecko_rsp_le_gap_set_advertise_report_scan_request_id         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x10030000)
#define gecko_rsp_le_gap_set_advertise_phy_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x11030000)
#define gecko_rsp_le_gap_set_advertise_configuration_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x12030000)
#define gecko_rsp_le_gap_clear_advertise_configuration_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x13030000)
#define gecko_rsp_le_gap_start_advertising_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x14030000)
#define gecko_rsp_le_gap_stop_advertising_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x15030000)
#define gecko_rsp_le_gap_set_discovery_timing_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x16030000)
#define gecko_rsp_le_gap_set_discovery_type_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x17030000)
#define gecko_rsp_le_gap_start_discovery_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x18030000)
#define gecko_rsp_le_gap_set_data_channel_classification_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x19030000)
#define gecko_rsp_le_gap_connect_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x1a030000)
#define gecko_rsp_le_gap_set_advertise_tx_power_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x1b030000)
#define gecko_rsp_le_gap_set_discovery_extended_scan_response_id      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x1c030000)
#define gecko_rsp_le_gap_start_periodic_advertising_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x1d030000)
#define gecko_rsp_le_gap_stop_periodic_advertising_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x1f030000)
#define gecko_rsp_le_gap_set_long_advertising_data_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x20030000)
#define gecko_rsp_le_gap_enable_whitelisting_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x21030000)
#define gecko_rsp_le_gap_set_conn_timing_parameters_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x22030000)
#define gecko_rsp_le_gap_set_advertise_random_address_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x25030000)
#define gecko_rsp_le_gap_clear_advertise_random_address_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x26030000)
#define gecko_rsp_sync_open_id                                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00420000)
#define gecko_rsp_sync_close_id                                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01420000)
#define gecko_rsp_le_connection_set_parameters_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00080000)
#define gecko_rsp_le_connection_get_rssi_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01080000)
#define gecko_rsp_le_connection_disable_slave_latency_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02080000)
#define gecko_rsp_le_connection_set_phy_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03080000)
#define gecko_rsp_le_connection_close_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04080000)
#define gecko_rsp_le_connection_set_timing_parameters_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05080000)
#define gecko_rsp_le_connection_read_channel_map_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x06080000)
#define gecko_rsp_le_connection_set_preferred_phy_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x07080000)
#define gecko_rsp_gatt_set_max_mtu_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00090000)
#define gecko_rsp_gatt_discover_primary_services_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01090000)
#define gecko_rsp_gatt_discover_primary_services_by_uuid_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02090000)
#define gecko_rsp_gatt_discover_characteristics_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03090000)
#define gecko_rsp_gatt_discover_characteristics_by_uuid_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04090000)
#define gecko_rsp_gatt_set_characteristic_notification_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05090000)
#define gecko_rsp_gatt_discover_descriptors_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x06090000)
#define gecko_rsp_gatt_read_characteristic_value_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x07090000)
#define gecko_rsp_gatt_read_characteristic_value_by_uuid_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x08090000)
#define gecko_rsp_gatt_write_characteristic_value_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x09090000)
#define gecko_rsp_gatt_write_characteristic_value_without_response_id  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0a090000)
#define gecko_rsp_gatt_prepare_characteristic_value_write_id          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0b090000)
#define gecko_rsp_gatt_execute_characteristic_value_write_id          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0c090000)
#define gecko_rsp_gatt_send_characteristic_confirmation_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0d090000)
#define gecko_rsp_gatt_read_descriptor_value_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0e090000)
#define gecko_rsp_gatt_write_descriptor_value_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0f090000)
#define gecko_rsp_gatt_find_included_services_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x10090000)
#define gecko_rsp_gatt_read_multiple_characteristic_values_id         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x11090000)
#define gecko_rsp_gatt_read_characteristic_value_from_offset_id       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x12090000)
#define gecko_rsp_gatt_prepare_characteristic_value_reliable_write_id  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x13090000)
#define gecko_rsp_gatt_server_read_attribute_value_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x000a0000)
#define gecko_rsp_gatt_server_read_attribute_type_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x010a0000)
#define gecko_rsp_gatt_server_write_attribute_value_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x020a0000)
#define gecko_rsp_gatt_server_send_user_read_response_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x030a0000)
#define gecko_rsp_gatt_server_send_user_write_response_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x040a0000)
#define gecko_rsp_gatt_server_send_characteristic_notification_id     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x050a0000)
#define gecko_rsp_gatt_server_find_attribute_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x060a0000)
#define gecko_rsp_gatt_server_set_capabilities_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x080a0000)
#define gecko_rsp_gatt_server_set_max_mtu_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0a0a0000)
#define gecko_rsp_gatt_server_get_mtu_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0b0a0000)
#define gecko_rsp_gatt_server_enable_capabilities_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0c0a0000)
#define gecko_rsp_gatt_server_disable_capabilities_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0d0a0000)
#define gecko_rsp_gatt_server_get_enabled_capabilities_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0e0a0000)
#define gecko_rsp_hardware_set_soft_timer_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x000c0000)
#define gecko_rsp_hardware_get_time_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0b0c0000)
#define gecko_rsp_hardware_set_lazy_soft_timer_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0c0c0000)
#define gecko_rsp_flash_ps_erase_all_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x010d0000)
#define gecko_rsp_flash_ps_save_id                                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x020d0000)
#define gecko_rsp_flash_ps_load_id                                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x030d0000)
#define gecko_rsp_flash_ps_erase_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x040d0000)
#define gecko_rsp_test_dtm_tx_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x000e0000)
#define gecko_rsp_test_dtm_rx_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x010e0000)
#define gecko_rsp_test_dtm_end_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x020e0000)
#define gecko_rsp_sm_set_bondable_mode_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x000f0000)
#define gecko_rsp_sm_configure_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x010f0000)
#define gecko_rsp_sm_store_bonding_configuration_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x020f0000)
#define gecko_rsp_sm_increase_security_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x040f0000)
#define gecko_rsp_sm_delete_bonding_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x060f0000)
#define gecko_rsp_sm_delete_bondings_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x070f0000)
#define gecko_rsp_sm_enter_passkey_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x080f0000)
#define gecko_rsp_sm_passkey_confirm_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x090f0000)
#define gecko_rsp_sm_set_oob_data_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0a0f0000)
#define gecko_rsp_sm_list_all_bondings_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0b0f0000)
#define gecko_rsp_sm_bonding_confirm_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0e0f0000)
#define gecko_rsp_sm_set_debug_mode_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0f0f0000)
#define gecko_rsp_sm_set_passkey_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x100f0000)
#define gecko_rsp_sm_use_sc_oob_id                                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x110f0000)
#define gecko_rsp_sm_set_sc_remote_oob_data_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x120f0000)
#define gecko_rsp_sm_add_to_whitelist_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x130f0000)
#define gecko_rsp_sm_set_minimum_key_size_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x140f0000)
#define gecko_rsp_homekit_configure_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00130000)
#define gecko_rsp_homekit_advertise_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01130000)
#define gecko_rsp_homekit_delete_pairings_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02130000)
#define gecko_rsp_homekit_check_authcp_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03130000)
#define gecko_rsp_homekit_get_pairing_id_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04130000)
#define gecko_rsp_homekit_send_write_response_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05130000)
#define gecko_rsp_homekit_send_read_response_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x06130000)
#define gecko_rsp_homekit_gsn_action_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x07130000)
#define gecko_rsp_homekit_event_notification_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x08130000)
#define gecko_rsp_homekit_broadcast_action_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x09130000)
#define gecko_rsp_homekit_configure_product_data_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x0a130000)
#define gecko_rsp_coex_set_options_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00200000)
#define gecko_rsp_coex_get_counters_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01200000)
#define gecko_rsp_coex_set_parameters_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02200000)
#define gecko_rsp_coex_set_directional_priority_pulse_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03200000)
#define gecko_rsp_l2cap_coc_send_connection_request_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01430000)
#define gecko_rsp_l2cap_coc_send_connection_response_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02430000)
#define gecko_rsp_l2cap_coc_send_le_flow_control_credit_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03430000)
#define gecko_rsp_l2cap_coc_send_disconnection_request_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04430000)
#define gecko_rsp_l2cap_coc_send_data_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05430000)
#define gecko_rsp_cte_transmitter_enable_connection_cte_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00440000)
#define gecko_rsp_cte_transmitter_disable_connection_cte_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01440000)
#define gecko_rsp_cte_transmitter_enable_connectionless_cte_id        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02440000)
#define gecko_rsp_cte_transmitter_disable_connectionless_cte_id       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03440000)
#define gecko_rsp_cte_transmitter_set_dtm_parameters_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04440000)
#define gecko_rsp_cte_transmitter_clear_dtm_parameters_id             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05440000)
#define gecko_rsp_cte_transmitter_enable_silabs_cte_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x06440000)
#define gecko_rsp_cte_transmitter_disable_silabs_cte_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x07440000)
#define gecko_rsp_cte_receiver_configure_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00450000)
#define gecko_rsp_cte_receiver_enable_connection_cte_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x01450000)
#define gecko_rsp_cte_receiver_disable_connection_cte_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x02450000)
#define gecko_rsp_cte_receiver_enable_connectionless_cte_id           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x03450000)
#define gecko_rsp_cte_receiver_disable_connectionless_cte_id          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x04450000)
#define gecko_rsp_cte_receiver_set_dtm_parameters_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x05450000)
#define gecko_rsp_cte_receiver_clear_dtm_parameters_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x06450000)
#define gecko_rsp_cte_receiver_enable_silabs_cte_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x07450000)
#define gecko_rsp_cte_receiver_disable_silabs_cte_id                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x08450000)
#define gecko_rsp_user_message_to_target_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_rsp|0x00ff0000)

#define gecko_evt_dfu_boot_id                                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00000000)
#define gecko_evt_dfu_boot_failure_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01000000)
#define gecko_evt_system_boot_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00010000)
#define gecko_evt_system_external_signal_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x03010000)
#define gecko_evt_system_awake_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x04010000)
#define gecko_evt_system_hardware_error_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x05010000)
#define gecko_evt_system_error_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x06010000)
#define gecko_evt_le_gap_scan_response_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00030000)
#define gecko_evt_le_gap_adv_timeout_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01030000)
#define gecko_evt_le_gap_scan_request_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02030000)
#define gecko_evt_le_gap_extended_scan_response_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x04030000)
#define gecko_evt_le_gap_periodic_advertising_status_id               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x05030000)
#define gecko_evt_sync_opened_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00420000)
#define gecko_evt_sync_closed_id                                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01420000)
#define gecko_evt_sync_data_id                                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02420000)
#define gecko_evt_le_connection_opened_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00080000)
#define gecko_evt_le_connection_closed_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01080000)
#define gecko_evt_le_connection_parameters_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02080000)
#define gecko_evt_le_connection_rssi_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x03080000)
#define gecko_evt_le_connection_phy_status_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x04080000)
#define gecko_evt_gatt_mtu_exchanged_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00090000)
#define gecko_evt_gatt_service_id                                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01090000)
#define gecko_evt_gatt_characteristic_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02090000)
#define gecko_evt_gatt_descriptor_id                                  (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x03090000)
#define gecko_evt_gatt_characteristic_value_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x04090000)
#define gecko_evt_gatt_descriptor_value_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x05090000)
#define gecko_evt_gatt_procedure_completed_id                         (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x06090000)
#define gecko_evt_gatt_server_attribute_value_id                      (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x000a0000)
#define gecko_evt_gatt_server_user_read_request_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x010a0000)
#define gecko_evt_gatt_server_user_write_request_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x020a0000)
#define gecko_evt_gatt_server_characteristic_status_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x030a0000)
#define gecko_evt_gatt_server_execute_write_completed_id              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x040a0000)
#define gecko_evt_hardware_soft_timer_id                              (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x000c0000)
#define gecko_evt_test_dtm_completed_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x000e0000)
#define gecko_evt_sm_passkey_display_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x000f0000)
#define gecko_evt_sm_passkey_request_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x010f0000)
#define gecko_evt_sm_confirm_passkey_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x020f0000)
#define gecko_evt_sm_bonded_id                                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x030f0000)
#define gecko_evt_sm_bonding_failed_id                                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x040f0000)
#define gecko_evt_sm_list_bonding_entry_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x050f0000)
#define gecko_evt_sm_list_all_bondings_complete_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x060f0000)
#define gecko_evt_sm_confirm_bonding_id                               (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x090f0000)
#define gecko_evt_homekit_setupcode_display_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00130000)
#define gecko_evt_homekit_paired_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01130000)
#define gecko_evt_homekit_pair_verified_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02130000)
#define gecko_evt_homekit_connection_opened_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x03130000)
#define gecko_evt_homekit_connection_closed_id                        (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x04130000)
#define gecko_evt_homekit_identify_id                                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x05130000)
#define gecko_evt_homekit_write_request_id                            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x06130000)
#define gecko_evt_homekit_read_request_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x07130000)
#define gecko_evt_homekit_disconnection_required_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x08130000)
#define gecko_evt_homekit_pairing_removed_id                          (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x09130000)
#define gecko_evt_homekit_setuppayload_display_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x0a130000)
#define gecko_evt_l2cap_coc_connection_request_id                     (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01430000)
#define gecko_evt_l2cap_coc_connection_response_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02430000)
#define gecko_evt_l2cap_coc_le_flow_control_credit_id                 (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x03430000)
#define gecko_evt_l2cap_coc_channel_disconnected_id                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x04430000)
#define gecko_evt_l2cap_coc_data_id                                   (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x05430000)
#define gecko_evt_l2cap_command_rejected_id                           (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x06430000)
#define gecko_evt_cte_receiver_connection_iq_report_id                (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00450000)
#define gecko_evt_cte_receiver_connectionless_iq_report_id            (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x01450000)
#define gecko_evt_cte_receiver_dtm_iq_report_id                       (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x02450000)
#define gecko_evt_cte_receiver_silabs_iq_report_id                    (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x03450000)
#define gecko_evt_user_message_to_host_id                             (((uint32)gecko_dev_type_gecko)|gecko_msg_type_evt|0x00ff0000)


PACKSTRUCT( struct gecko_msg_dfu_reset_cmd_t
{
    uint8               dfu;
});
PACKSTRUCT( struct gecko_msg_dfu_flash_set_address_cmd_t
{
    uint32              address;
});
PACKSTRUCT( struct gecko_msg_dfu_flash_set_address_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_dfu_flash_upload_cmd_t
{
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_dfu_flash_upload_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_dfu_flash_upload_finish_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_dfu_boot_evt_t
{
    uint32              version;
});
PACKSTRUCT( struct gecko_msg_dfu_boot_failure_evt_t
{
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_system_hello_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_reset_cmd_t
{
    uint8               dfu;
});
PACKSTRUCT( struct gecko_msg_system_get_bt_address_rsp_t
{
    bd_addr             address;
});
PACKSTRUCT( struct gecko_msg_system_set_bt_address_cmd_t
{
    bd_addr             address;
});
PACKSTRUCT( struct gecko_msg_system_set_bt_address_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_set_tx_power_cmd_t
{
    int16               power;
});
PACKSTRUCT( struct gecko_msg_system_set_tx_power_rsp_t
{
    int16               set_power;
});
PACKSTRUCT( struct gecko_msg_system_get_random_data_cmd_t
{
    uint8               length;
});
PACKSTRUCT( struct gecko_msg_system_get_random_data_rsp_t
{
    uint16              result;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_system_halt_cmd_t
{
    uint8               halt;
});
PACKSTRUCT( struct gecko_msg_system_halt_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_set_device_name_cmd_t
{
    uint8               type;
    uint8array          name;
});
PACKSTRUCT( struct gecko_msg_system_set_device_name_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_linklayer_configure_cmd_t
{
    uint8               key;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_system_linklayer_configure_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_get_counters_cmd_t
{
    uint8               reset;
});
PACKSTRUCT( struct gecko_msg_system_get_counters_rsp_t
{
    uint16              result;
    uint16              tx_packets;
    uint16              rx_packets;
    uint16              crc_errors;
    uint16              failures;
});
PACKSTRUCT( struct gecko_msg_system_data_buffer_write_cmd_t
{
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_system_data_buffer_write_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_set_identity_address_cmd_t
{
    bd_addr             address;
    uint8               type;
});
PACKSTRUCT( struct gecko_msg_system_set_identity_address_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_data_buffer_clear_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_system_boot_evt_t
{
    uint16              major;
    uint16              minor;
    uint16              patch;
    uint16              build;
    uint32              bootloader;
    uint16              hw;
    uint32              hash;
});
PACKSTRUCT( struct gecko_msg_system_external_signal_evt_t
{
    uint32              extsignals;
});
PACKSTRUCT( struct gecko_msg_system_hardware_error_evt_t
{
    uint16              status;
});
PACKSTRUCT( struct gecko_msg_system_error_evt_t
{
    uint16              reason;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_le_gap_open_cmd_t
{
    bd_addr             address;
    uint8               address_type;
});
PACKSTRUCT( struct gecko_msg_le_gap_open_rsp_t
{
    uint16              result;
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_mode_cmd_t
{
    uint8               discover;
    uint8               connect;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_mode_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_discover_cmd_t
{
    uint8               mode;
});
PACKSTRUCT( struct gecko_msg_le_gap_discover_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_end_procedure_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_adv_parameters_cmd_t
{
    uint16              interval_min;
    uint16              interval_max;
    uint8               channel_map;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_adv_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_conn_parameters_cmd_t
{
    uint16              min_interval;
    uint16              max_interval;
    uint16              latency;
    uint16              timeout;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_conn_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_scan_parameters_cmd_t
{
    uint16              scan_interval;
    uint16              scan_window;
    uint8               active;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_scan_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_adv_data_cmd_t
{
    uint8               scan_rsp;
    uint8array          adv_data;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_adv_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_adv_timeout_cmd_t
{
    uint8               maxevents;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_adv_timeout_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_conn_phy_cmd_t
{
    uint8               preferred_phy;
    uint8               accepted_phy;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_conn_phy_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_bt5_set_mode_cmd_t
{
    uint8               handle;
    uint8               discover;
    uint8               connect;
    uint16              maxevents;
    uint8               address_type;
});
PACKSTRUCT( struct gecko_msg_le_gap_bt5_set_mode_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_bt5_set_adv_parameters_cmd_t
{
    uint8               handle;
    uint16              interval_min;
    uint16              interval_max;
    uint8               channel_map;
    uint8               report_scan;
});
PACKSTRUCT( struct gecko_msg_le_gap_bt5_set_adv_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_bt5_set_adv_data_cmd_t
{
    uint8               handle;
    uint8               scan_rsp;
    uint8array          adv_data;
});
PACKSTRUCT( struct gecko_msg_le_gap_bt5_set_adv_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_privacy_mode_cmd_t
{
    uint8               privacy;
    uint8               interval;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_privacy_mode_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_timing_cmd_t
{
    uint8               handle;
    uint32              interval_min;
    uint32              interval_max;
    uint16              duration;
    uint8               maxevents;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_timing_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_channel_map_cmd_t
{
    uint8               handle;
    uint8               channel_map;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_channel_map_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_report_scan_request_cmd_t
{
    uint8               handle;
    uint8               report_scan_req;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_report_scan_request_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_phy_cmd_t
{
    uint8               handle;
    uint8               primary_phy;
    uint8               secondary_phy;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_phy_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_configuration_cmd_t
{
    uint8               handle;
    uint32              configurations;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_configuration_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_clear_advertise_configuration_cmd_t
{
    uint8               handle;
    uint32              configurations;
});
PACKSTRUCT( struct gecko_msg_le_gap_clear_advertise_configuration_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_start_advertising_cmd_t
{
    uint8               handle;
    uint8               discover;
    uint8               connect;
});
PACKSTRUCT( struct gecko_msg_le_gap_start_advertising_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_stop_advertising_cmd_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_le_gap_stop_advertising_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_discovery_timing_cmd_t
{
    uint8               phys;
    uint16              scan_interval;
    uint16              scan_window;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_discovery_timing_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_discovery_type_cmd_t
{
    uint8               phys;
    uint8               scan_type;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_discovery_type_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_start_discovery_cmd_t
{
    uint8               scanning_phy;
    uint8               mode;
});
PACKSTRUCT( struct gecko_msg_le_gap_start_discovery_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_data_channel_classification_cmd_t
{
    uint8array          channel_map;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_data_channel_classification_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_connect_cmd_t
{
    bd_addr             address;
    uint8               address_type;
    uint8               initiating_phy;
});
PACKSTRUCT( struct gecko_msg_le_gap_connect_rsp_t
{
    uint16              result;
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_tx_power_cmd_t
{
    uint8               handle;
    int16               power;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_tx_power_rsp_t
{
    uint16              result;
    int16               set_power;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_discovery_extended_scan_response_cmd_t
{
    uint8               enable;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_discovery_extended_scan_response_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_start_periodic_advertising_cmd_t
{
    uint8               handle;
    uint16              interval_min;
    uint16              interval_max;
    uint32              flags;
});
PACKSTRUCT( struct gecko_msg_le_gap_start_periodic_advertising_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_stop_periodic_advertising_cmd_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_le_gap_stop_periodic_advertising_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_long_advertising_data_cmd_t
{
    uint8               handle;
    uint8               packet_type;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_long_advertising_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_enable_whitelisting_cmd_t
{
    uint8               enable;
});
PACKSTRUCT( struct gecko_msg_le_gap_enable_whitelisting_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_conn_timing_parameters_cmd_t
{
    uint16              min_interval;
    uint16              max_interval;
    uint16              latency;
    uint16              timeout;
    uint16              min_ce_length;
    uint16              max_ce_length;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_conn_timing_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_random_address_cmd_t
{
    uint8               handle;
    uint8               addr_type;
    bd_addr             address;
});
PACKSTRUCT( struct gecko_msg_le_gap_set_advertise_random_address_rsp_t
{
    uint16              result;
    bd_addr             address_out;
});
PACKSTRUCT( struct gecko_msg_le_gap_clear_advertise_random_address_cmd_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_le_gap_clear_advertise_random_address_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_gap_scan_response_evt_t
{
    int8                rssi;
    uint8               packet_type;
    bd_addr             address;
    uint8               address_type;
    uint8               bonding;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_le_gap_adv_timeout_evt_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_le_gap_scan_request_evt_t
{
    uint8               handle;
    bd_addr             address;
    uint8               address_type;
    uint8               bonding;
});
PACKSTRUCT( struct gecko_msg_le_gap_extended_scan_response_evt_t
{
    uint8               packet_type;
    bd_addr             address;
    uint8               address_type;
    uint8               bonding;
    uint8               primary_phy;
    uint8               secondary_phy;
    uint8               adv_sid;
    int8                tx_power;
    int8                rssi;
    uint8               channel;
    uint16              periodic_interval;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_le_gap_periodic_advertising_status_evt_t
{
    uint8               sid;
    uint32              status;
});
PACKSTRUCT( struct gecko_msg_sync_open_cmd_t
{
    uint8               adv_sid;
    uint16              skip;
    uint16              timeout;
    bd_addr             address;
    uint8               address_type;
});
PACKSTRUCT( struct gecko_msg_sync_open_rsp_t
{
    uint16              result;
    uint8               sync;
});
PACKSTRUCT( struct gecko_msg_sync_close_cmd_t
{
    uint8               sync;
});
PACKSTRUCT( struct gecko_msg_sync_close_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sync_opened_evt_t
{
    uint8               sync;
    uint8               adv_sid;
    bd_addr             address;
    uint8               address_type;
    uint8               adv_phy;
    uint16              adv_interval;
    uint16              clock_accuracy;
});
PACKSTRUCT( struct gecko_msg_sync_closed_evt_t
{
    uint16              reason;
    uint8               sync;
});
PACKSTRUCT( struct gecko_msg_sync_data_evt_t
{
    uint8               sync;
    int8                tx_power;
    int8                rssi;
    uint8               data_status;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_parameters_cmd_t
{
    uint8               connection;
    uint16              min_interval;
    uint16              max_interval;
    uint16              latency;
    uint16              timeout;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_get_rssi_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_le_connection_get_rssi_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_disable_slave_latency_cmd_t
{
    uint8               connection;
    uint8               disable;
});
PACKSTRUCT( struct gecko_msg_le_connection_disable_slave_latency_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_phy_cmd_t
{
    uint8               connection;
    uint8               phy;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_phy_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_close_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_le_connection_close_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_timing_parameters_cmd_t
{
    uint8               connection;
    uint16              min_interval;
    uint16              max_interval;
    uint16              latency;
    uint16              timeout;
    uint16              min_ce_length;
    uint16              max_ce_length;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_timing_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_read_channel_map_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_le_connection_read_channel_map_rsp_t
{
    uint16              result;
    uint8array          channel_map;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_preferred_phy_cmd_t
{
    uint8               connection;
    uint8               preferred_phy;
    uint8               accepted_phy;
});
PACKSTRUCT( struct gecko_msg_le_connection_set_preferred_phy_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_le_connection_opened_evt_t
{
    bd_addr             address;
    uint8               address_type;
    uint8               master;
    uint8               connection;
    uint8               bonding;
    uint8               advertiser;
});
PACKSTRUCT( struct gecko_msg_le_connection_closed_evt_t
{
    uint16              reason;
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_le_connection_parameters_evt_t
{
    uint8               connection;
    uint16              interval;
    uint16              latency;
    uint16              timeout;
    uint8               security_mode;
    uint16              txsize;
});
PACKSTRUCT( struct gecko_msg_le_connection_rssi_evt_t
{
    uint8               connection;
    uint8               status;
    int8                rssi;
});
PACKSTRUCT( struct gecko_msg_le_connection_phy_status_evt_t
{
    uint8               connection;
    uint8               phy;
});
PACKSTRUCT( struct gecko_msg_gatt_set_max_mtu_cmd_t
{
    uint16              max_mtu;
});
PACKSTRUCT( struct gecko_msg_gatt_set_max_mtu_rsp_t
{
    uint16              result;
    uint16              max_mtu;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_primary_services_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_primary_services_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_primary_services_by_uuid_cmd_t
{
    uint8               connection;
    uint8array          uuid;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_primary_services_by_uuid_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_characteristics_cmd_t
{
    uint8               connection;
    uint32              service;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_characteristics_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_characteristics_by_uuid_cmd_t
{
    uint8               connection;
    uint32              service;
    uint8array          uuid;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_characteristics_by_uuid_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_set_characteristic_notification_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               flags;
});
PACKSTRUCT( struct gecko_msg_gatt_set_characteristic_notification_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_descriptors_cmd_t
{
    uint8               connection;
    uint16              characteristic;
});
PACKSTRUCT( struct gecko_msg_gatt_discover_descriptors_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_read_characteristic_value_cmd_t
{
    uint8               connection;
    uint16              characteristic;
});
PACKSTRUCT( struct gecko_msg_gatt_read_characteristic_value_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_read_characteristic_value_by_uuid_cmd_t
{
    uint8               connection;
    uint32              service;
    uint8array          uuid;
});
PACKSTRUCT( struct gecko_msg_gatt_read_characteristic_value_by_uuid_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_write_characteristic_value_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_write_characteristic_value_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_write_characteristic_value_without_response_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_write_characteristic_value_without_response_rsp_t
{
    uint16              result;
    uint16              sent_len;
});
PACKSTRUCT( struct gecko_msg_gatt_prepare_characteristic_value_write_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_prepare_characteristic_value_write_rsp_t
{
    uint16              result;
    uint16              sent_len;
});
PACKSTRUCT( struct gecko_msg_gatt_execute_characteristic_value_write_cmd_t
{
    uint8               connection;
    uint8               flags;
});
PACKSTRUCT( struct gecko_msg_gatt_execute_characteristic_value_write_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_send_characteristic_confirmation_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_gatt_send_characteristic_confirmation_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_read_descriptor_value_cmd_t
{
    uint8               connection;
    uint16              descriptor;
});
PACKSTRUCT( struct gecko_msg_gatt_read_descriptor_value_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_write_descriptor_value_cmd_t
{
    uint8               connection;
    uint16              descriptor;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_write_descriptor_value_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_find_included_services_cmd_t
{
    uint8               connection;
    uint32              service;
});
PACKSTRUCT( struct gecko_msg_gatt_find_included_services_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_read_multiple_characteristic_values_cmd_t
{
    uint8               connection;
    uint8array          characteristic_list;
});
PACKSTRUCT( struct gecko_msg_gatt_read_multiple_characteristic_values_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_read_characteristic_value_from_offset_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint16              offset;
    uint16              maxlen;
});
PACKSTRUCT( struct gecko_msg_gatt_read_characteristic_value_from_offset_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_prepare_characteristic_value_reliable_write_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_prepare_characteristic_value_reliable_write_rsp_t
{
    uint16              result;
    uint16              sent_len;
});
PACKSTRUCT( struct gecko_msg_gatt_mtu_exchanged_evt_t
{
    uint8               connection;
    uint16              mtu;
});
PACKSTRUCT( struct gecko_msg_gatt_service_evt_t
{
    uint8               connection;
    uint32              service;
    uint8array          uuid;
});
PACKSTRUCT( struct gecko_msg_gatt_characteristic_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               properties;
    uint8array          uuid;
});
PACKSTRUCT( struct gecko_msg_gatt_descriptor_evt_t
{
    uint8               connection;
    uint16              descriptor;
    uint8array          uuid;
});
PACKSTRUCT( struct gecko_msg_gatt_characteristic_value_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               att_opcode;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_descriptor_value_evt_t
{
    uint8               connection;
    uint16              descriptor;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_procedure_completed_evt_t
{
    uint8               connection;
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_server_read_attribute_value_cmd_t
{
    uint16              attribute;
    uint16              offset;
});
PACKSTRUCT( struct gecko_msg_gatt_server_read_attribute_value_rsp_t
{
    uint16              result;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_server_read_attribute_type_cmd_t
{
    uint16              attribute;
});
PACKSTRUCT( struct gecko_msg_gatt_server_read_attribute_type_rsp_t
{
    uint16              result;
    uint8array          type;
});
PACKSTRUCT( struct gecko_msg_gatt_server_write_attribute_value_cmd_t
{
    uint16              attribute;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_server_write_attribute_value_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_server_send_user_read_response_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               att_errorcode;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_server_send_user_read_response_rsp_t
{
    uint16              result;
    uint16              sent_len;
});
PACKSTRUCT( struct gecko_msg_gatt_server_send_user_write_response_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               att_errorcode;
});
PACKSTRUCT( struct gecko_msg_gatt_server_send_user_write_response_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_server_send_characteristic_notification_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_server_send_characteristic_notification_rsp_t
{
    uint16              result;
    uint16              sent_len;
});
PACKSTRUCT( struct gecko_msg_gatt_server_find_attribute_cmd_t
{
    uint16              start;
    uint8array          type;
});
PACKSTRUCT( struct gecko_msg_gatt_server_find_attribute_rsp_t
{
    uint16              result;
    uint16              attribute;
});
PACKSTRUCT( struct gecko_msg_gatt_server_set_capabilities_cmd_t
{
    uint32              caps;
    uint32              reserved;
});
PACKSTRUCT( struct gecko_msg_gatt_server_set_capabilities_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_server_set_max_mtu_cmd_t
{
    uint16              max_mtu;
});
PACKSTRUCT( struct gecko_msg_gatt_server_set_max_mtu_rsp_t
{
    uint16              result;
    uint16              max_mtu;
});
PACKSTRUCT( struct gecko_msg_gatt_server_get_mtu_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_gatt_server_get_mtu_rsp_t
{
    uint16              result;
    uint16              mtu;
});
PACKSTRUCT( struct gecko_msg_gatt_server_enable_capabilities_cmd_t
{
    uint32              caps;
});
PACKSTRUCT( struct gecko_msg_gatt_server_enable_capabilities_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_server_disable_capabilities_cmd_t
{
    uint32              caps;
});
PACKSTRUCT( struct gecko_msg_gatt_server_disable_capabilities_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_gatt_server_get_enabled_capabilities_rsp_t
{
    uint16              result;
    uint32              caps;
});
PACKSTRUCT( struct gecko_msg_gatt_server_attribute_value_evt_t
{
    uint8               connection;
    uint16              attribute;
    uint8               att_opcode;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_server_user_read_request_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               att_opcode;
    uint16              offset;
});
PACKSTRUCT( struct gecko_msg_gatt_server_user_write_request_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               att_opcode;
    uint16              offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_gatt_server_characteristic_status_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               status_flags;
    uint16              client_config_flags;
});
PACKSTRUCT( struct gecko_msg_gatt_server_execute_write_completed_evt_t
{
    uint8               connection;
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_hardware_set_soft_timer_cmd_t
{
    uint32              time;
    uint8               handle;
    uint8               single_shot;
});
PACKSTRUCT( struct gecko_msg_hardware_set_soft_timer_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_hardware_get_time_rsp_t
{
    uint32              seconds;
    uint16              ticks;
});
PACKSTRUCT( struct gecko_msg_hardware_set_lazy_soft_timer_cmd_t
{
    uint32              time;
    uint32              slack;
    uint8               handle;
    uint8               single_shot;
});
PACKSTRUCT( struct gecko_msg_hardware_set_lazy_soft_timer_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_hardware_soft_timer_evt_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_flash_ps_erase_all_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_flash_ps_save_cmd_t
{
    uint16              key;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_flash_ps_save_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_flash_ps_load_cmd_t
{
    uint16              key;
});
PACKSTRUCT( struct gecko_msg_flash_ps_load_rsp_t
{
    uint16              result;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_flash_ps_erase_cmd_t
{
    uint16              key;
});
PACKSTRUCT( struct gecko_msg_flash_ps_erase_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_test_dtm_tx_cmd_t
{
    uint8               packet_type;
    uint8               length;
    uint8               channel;
    uint8               phy;
});
PACKSTRUCT( struct gecko_msg_test_dtm_tx_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_test_dtm_rx_cmd_t
{
    uint8               channel;
    uint8               phy;
});
PACKSTRUCT( struct gecko_msg_test_dtm_rx_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_test_dtm_end_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_test_dtm_completed_evt_t
{
    uint16              result;
    uint16              number_of_packets;
});
PACKSTRUCT( struct gecko_msg_sm_set_bondable_mode_cmd_t
{
    uint8               bondable;
});
PACKSTRUCT( struct gecko_msg_sm_set_bondable_mode_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_configure_cmd_t
{
    uint8               flags;
    uint8               io_capabilities;
});
PACKSTRUCT( struct gecko_msg_sm_configure_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_store_bonding_configuration_cmd_t
{
    uint8               max_bonding_count;
    uint8               policy_flags;
});
PACKSTRUCT( struct gecko_msg_sm_store_bonding_configuration_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_increase_security_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_sm_increase_security_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_delete_bonding_cmd_t
{
    uint8               bonding;
});
PACKSTRUCT( struct gecko_msg_sm_delete_bonding_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_delete_bondings_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_enter_passkey_cmd_t
{
    uint8               connection;
    int32               passkey;
});
PACKSTRUCT( struct gecko_msg_sm_enter_passkey_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_passkey_confirm_cmd_t
{
    uint8               connection;
    uint8               confirm;
});
PACKSTRUCT( struct gecko_msg_sm_passkey_confirm_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_set_oob_data_cmd_t
{
    uint8array          oob_data;
});
PACKSTRUCT( struct gecko_msg_sm_set_oob_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_list_all_bondings_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_bonding_confirm_cmd_t
{
    uint8               connection;
    uint8               confirm;
});
PACKSTRUCT( struct gecko_msg_sm_bonding_confirm_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_set_debug_mode_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_set_passkey_cmd_t
{
    int32               passkey;
});
PACKSTRUCT( struct gecko_msg_sm_set_passkey_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_use_sc_oob_cmd_t
{
    uint8               enable;
});
PACKSTRUCT( struct gecko_msg_sm_use_sc_oob_rsp_t
{
    uint16              result;
    uint8array          oob_data;
});
PACKSTRUCT( struct gecko_msg_sm_set_sc_remote_oob_data_cmd_t
{
    uint8array          oob_data;
});
PACKSTRUCT( struct gecko_msg_sm_set_sc_remote_oob_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_add_to_whitelist_cmd_t
{
    bd_addr             address;
    uint8               address_type;
});
PACKSTRUCT( struct gecko_msg_sm_add_to_whitelist_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_set_minimum_key_size_cmd_t
{
    uint8               minimum_key_size;
});
PACKSTRUCT( struct gecko_msg_sm_set_minimum_key_size_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_sm_passkey_display_evt_t
{
    uint8               connection;
    uint32              passkey;
});
PACKSTRUCT( struct gecko_msg_sm_passkey_request_evt_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_sm_confirm_passkey_evt_t
{
    uint8               connection;
    uint32              passkey;
});
PACKSTRUCT( struct gecko_msg_sm_bonded_evt_t
{
    uint8               connection;
    uint8               bonding;
});
PACKSTRUCT( struct gecko_msg_sm_bonding_failed_evt_t
{
    uint8               connection;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_sm_list_bonding_entry_evt_t
{
    uint8               bonding;
    bd_addr             address;
    uint8               address_type;
});
PACKSTRUCT( struct gecko_msg_sm_confirm_bonding_evt_t
{
    uint8               connection;
    int8                bonding_handle;
});
PACKSTRUCT( struct gecko_msg_homekit_configure_cmd_t
{
    uint8               i2c_address;
    uint8               support_display;
    uint8               hap_attribute_features;
    uint16              category;
    uint8               configuration_number;
    uint16              fast_advert_interval;
    uint16              fast_advert_timeout;
    uint32              flag;
    uint16              broadcast_advert_timeout;
    uint8array          model_name;
});
PACKSTRUCT( struct gecko_msg_homekit_configure_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_advertise_cmd_t
{
    uint8               enable;
    uint16              interval_min;
    uint16              interval_max;
    uint8               channel_map;
});
PACKSTRUCT( struct gecko_msg_homekit_advertise_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_delete_pairings_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_check_authcp_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_get_pairing_id_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_homekit_get_pairing_id_rsp_t
{
    uint16              result;
    uint8array          pairing_id;
});
PACKSTRUCT( struct gecko_msg_homekit_send_write_response_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               status_code;
});
PACKSTRUCT( struct gecko_msg_homekit_send_write_response_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_send_read_response_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               status_code;
    uint16              attribute_size;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_homekit_send_read_response_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_gsn_action_cmd_t
{
    uint8               action;
});
PACKSTRUCT( struct gecko_msg_homekit_gsn_action_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_event_notification_cmd_t
{
    uint8               connection;
    uint16              characteristic;
    uint8               change_originator;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_homekit_event_notification_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_broadcast_action_cmd_t
{
    uint8               action;
    uint8array          params;
});
PACKSTRUCT( struct gecko_msg_homekit_broadcast_action_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_configure_product_data_cmd_t
{
    uint8array          product_data;
});
PACKSTRUCT( struct gecko_msg_homekit_configure_product_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_homekit_setupcode_display_evt_t
{
    uint8               connection;
    uint8array          setupcode;
});
PACKSTRUCT( struct gecko_msg_homekit_paired_evt_t
{
    uint8               connection;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_homekit_pair_verified_evt_t
{
    uint8               connection;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_homekit_connection_opened_evt_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_homekit_connection_closed_evt_t
{
    uint8               connection;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_homekit_identify_evt_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_homekit_write_request_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint16              chr_value_size;
    uint16              authorization_size;
    uint16              value_offset;
    uint8array          value;
});
PACKSTRUCT( struct gecko_msg_homekit_read_request_evt_t
{
    uint8               connection;
    uint16              characteristic;
    uint16              offset;
});
PACKSTRUCT( struct gecko_msg_homekit_disconnection_required_evt_t
{
    uint8               connection;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_homekit_pairing_removed_evt_t
{
    uint8               connection;
    uint16              remaining_pairings;
    uint8array          pairing_id;
});
PACKSTRUCT( struct gecko_msg_homekit_setuppayload_display_evt_t
{
    uint8               connection;
    uint8array          setuppayload;
});
PACKSTRUCT( struct gecko_msg_coex_set_options_cmd_t
{
    uint32              mask;
    uint32              options;
});
PACKSTRUCT( struct gecko_msg_coex_set_options_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_coex_get_counters_cmd_t
{
    uint8               reset;
});
PACKSTRUCT( struct gecko_msg_coex_get_counters_rsp_t
{
    uint16              result;
    uint8array          counters;
});
PACKSTRUCT( struct gecko_msg_coex_set_parameters_cmd_t
{
    uint8               priority;
    uint8               request;
    uint8               pwm_period;
    uint8               pwm_dutycycle;
});
PACKSTRUCT( struct gecko_msg_coex_set_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_coex_set_directional_priority_pulse_cmd_t
{
    uint8               pulse;
});
PACKSTRUCT( struct gecko_msg_coex_set_directional_priority_pulse_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_connection_request_cmd_t
{
    uint8               connection;
    uint16              le_psm;
    uint16              mtu;
    uint16              mps;
    uint16              initial_credit;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_connection_request_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_connection_response_cmd_t
{
    uint8               connection;
    uint16              cid;
    uint16              mtu;
    uint16              mps;
    uint16              initial_credit;
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_connection_response_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_le_flow_control_credit_cmd_t
{
    uint8               connection;
    uint16              cid;
    uint16              credits;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_le_flow_control_credit_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_disconnection_request_cmd_t
{
    uint8               connection;
    uint16              cid;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_disconnection_request_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_data_cmd_t
{
    uint8               connection;
    uint16              cid;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_send_data_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_connection_request_evt_t
{
    uint8               connection;
    uint16              le_psm;
    uint16              source_cid;
    uint16              mtu;
    uint16              mps;
    uint16              initial_credit;
    uint8               flags;
    uint8               encryption_key_size;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_connection_response_evt_t
{
    uint8               connection;
    uint16              destination_cid;
    uint16              mtu;
    uint16              mps;
    uint16              initial_credit;
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_le_flow_control_credit_evt_t
{
    uint8               connection;
    uint16              cid;
    uint16              credits;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_channel_disconnected_evt_t
{
    uint8               connection;
    uint16              cid;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_l2cap_coc_data_evt_t
{
    uint8               connection;
    uint16              cid;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_l2cap_command_rejected_evt_t
{
    uint8               connection;
    uint8               code;
    uint16              reason;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_enable_connection_cte_cmd_t
{
    uint8               connection;
    uint8               cte_types;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_enable_connection_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_disable_connection_cte_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_disable_connection_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_enable_connectionless_cte_cmd_t
{
    uint8               handle;
    uint8               cte_length;
    uint8               cte_type;
    uint8               cte_count;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_enable_connectionless_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_disable_connectionless_cte_cmd_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_disable_connectionless_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_set_dtm_parameters_cmd_t
{
    uint8               cte_length;
    uint8               cte_type;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_set_dtm_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_clear_dtm_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_enable_silabs_cte_cmd_t
{
    uint8               handle;
    uint8               cte_length;
    uint8               cte_type;
    uint8               cte_count;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_enable_silabs_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_disable_silabs_cte_cmd_t
{
    uint8               handle;
});
PACKSTRUCT( struct gecko_msg_cte_transmitter_disable_silabs_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_configure_cmd_t
{
    uint8               flags;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_configure_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_enable_connection_cte_cmd_t
{
    uint8               connection;
    uint16              interval;
    uint8               cte_length;
    uint8               cte_type;
    uint8               slot_durations;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_enable_connection_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_disable_connection_cte_cmd_t
{
    uint8               connection;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_disable_connection_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_enable_connectionless_cte_cmd_t
{
    uint8               sync;
    uint8               slot_durations;
    uint8               cte_count;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_enable_connectionless_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_disable_connectionless_cte_cmd_t
{
    uint8               sync;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_disable_connectionless_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_set_dtm_parameters_cmd_t
{
    uint8               cte_length;
    uint8               cte_type;
    uint8               slot_durations;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_set_dtm_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_clear_dtm_parameters_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_enable_silabs_cte_cmd_t
{
    uint8               slot_durations;
    uint8               cte_count;
    uint8array          switching_pattern;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_enable_silabs_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_disable_silabs_cte_rsp_t
{
    uint16              result;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_connection_iq_report_evt_t
{
    uint16              status;
    uint8               connection;
    uint8               phy;
    uint8               channel;
    int8                rssi;
    uint8               rssi_antenna_id;
    uint8               cte_type;
    uint8               slot_durations;
    uint16              event_counter;
    uint8array          samples;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_connectionless_iq_report_evt_t
{
    uint16              status;
    uint8               sync;
    uint8               channel;
    int8                rssi;
    uint8               rssi_antenna_id;
    uint8               cte_type;
    uint8               slot_durations;
    uint16              event_counter;
    uint8array          samples;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_dtm_iq_report_evt_t
{
    uint16              status;
    uint8               channel;
    int8                rssi;
    uint8               rssi_antenna_id;
    uint8               cte_type;
    uint8               slot_durations;
    uint16              event_counter;
    uint8array          samples;
});
PACKSTRUCT( struct gecko_msg_cte_receiver_silabs_iq_report_evt_t
{
    uint16              status;
    bd_addr             address;
    uint8               address_type;
    uint8               phy;
    uint8               channel;
    int8                rssi;
    uint8               rssi_antenna_id;
    uint8               cte_type;
    uint8               slot_durations;
    uint16              packet_counter;
    uint8array          samples;
});
PACKSTRUCT( struct gecko_msg_user_message_to_target_cmd_t
{
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_user_message_to_target_rsp_t
{
    uint16              result;
    uint8array          data;
});
PACKSTRUCT( struct gecko_msg_user_message_to_host_evt_t
{
    uint8array          data;
});


PACKSTRUCT( struct gecko_cmd_packet
{
    uint32   header;

union{
    uint8 handle;
    struct gecko_msg_dfu_reset_cmd_t                             cmd_dfu_reset;
    struct gecko_msg_dfu_flash_set_address_cmd_t                 cmd_dfu_flash_set_address;
    struct gecko_msg_dfu_flash_set_address_rsp_t                 rsp_dfu_flash_set_address;
    struct gecko_msg_dfu_flash_upload_cmd_t                      cmd_dfu_flash_upload;
    struct gecko_msg_dfu_flash_upload_rsp_t                      rsp_dfu_flash_upload;
    struct gecko_msg_dfu_flash_upload_finish_rsp_t               rsp_dfu_flash_upload_finish;
    struct gecko_msg_dfu_boot_evt_t                              evt_dfu_boot;
    struct gecko_msg_dfu_boot_failure_evt_t                      evt_dfu_boot_failure;
    struct gecko_msg_system_hello_rsp_t                          rsp_system_hello;
    struct gecko_msg_system_reset_cmd_t                          cmd_system_reset;
    struct gecko_msg_system_get_bt_address_rsp_t                 rsp_system_get_bt_address;
    struct gecko_msg_system_set_bt_address_cmd_t                 cmd_system_set_bt_address;
    struct gecko_msg_system_set_bt_address_rsp_t                 rsp_system_set_bt_address;
    struct gecko_msg_system_set_tx_power_cmd_t                   cmd_system_set_tx_power;
    struct gecko_msg_system_set_tx_power_rsp_t                   rsp_system_set_tx_power;
    struct gecko_msg_system_get_random_data_cmd_t                cmd_system_get_random_data;
    struct gecko_msg_system_get_random_data_rsp_t                rsp_system_get_random_data;
    struct gecko_msg_system_halt_cmd_t                           cmd_system_halt;
    struct gecko_msg_system_halt_rsp_t                           rsp_system_halt;
    struct gecko_msg_system_set_device_name_cmd_t                cmd_system_set_device_name;
    struct gecko_msg_system_set_device_name_rsp_t                rsp_system_set_device_name;
    struct gecko_msg_system_linklayer_configure_cmd_t            cmd_system_linklayer_configure;
    struct gecko_msg_system_linklayer_configure_rsp_t            rsp_system_linklayer_configure;
    struct gecko_msg_system_get_counters_cmd_t                   cmd_system_get_counters;
    struct gecko_msg_system_get_counters_rsp_t                   rsp_system_get_counters;
    struct gecko_msg_system_data_buffer_write_cmd_t              cmd_system_data_buffer_write;
    struct gecko_msg_system_data_buffer_write_rsp_t              rsp_system_data_buffer_write;
    struct gecko_msg_system_set_identity_address_cmd_t           cmd_system_set_identity_address;
    struct gecko_msg_system_set_identity_address_rsp_t           rsp_system_set_identity_address;
    struct gecko_msg_system_data_buffer_clear_rsp_t              rsp_system_data_buffer_clear;
    struct gecko_msg_system_boot_evt_t                           evt_system_boot;
    struct gecko_msg_system_external_signal_evt_t                evt_system_external_signal;
    struct gecko_msg_system_hardware_error_evt_t                 evt_system_hardware_error;
    struct gecko_msg_system_error_evt_t                          evt_system_error;
    struct gecko_msg_le_gap_open_cmd_t                           cmd_le_gap_open;
    struct gecko_msg_le_gap_open_rsp_t                           rsp_le_gap_open;
    struct gecko_msg_le_gap_set_mode_cmd_t                       cmd_le_gap_set_mode;
    struct gecko_msg_le_gap_set_mode_rsp_t                       rsp_le_gap_set_mode;
    struct gecko_msg_le_gap_discover_cmd_t                       cmd_le_gap_discover;
    struct gecko_msg_le_gap_discover_rsp_t                       rsp_le_gap_discover;
    struct gecko_msg_le_gap_end_procedure_rsp_t                  rsp_le_gap_end_procedure;
    struct gecko_msg_le_gap_set_adv_parameters_cmd_t             cmd_le_gap_set_adv_parameters;
    struct gecko_msg_le_gap_set_adv_parameters_rsp_t             rsp_le_gap_set_adv_parameters;
    struct gecko_msg_le_gap_set_conn_parameters_cmd_t            cmd_le_gap_set_conn_parameters;
    struct gecko_msg_le_gap_set_conn_parameters_rsp_t            rsp_le_gap_set_conn_parameters;
    struct gecko_msg_le_gap_set_scan_parameters_cmd_t            cmd_le_gap_set_scan_parameters;
    struct gecko_msg_le_gap_set_scan_parameters_rsp_t            rsp_le_gap_set_scan_parameters;
    struct gecko_msg_le_gap_set_adv_data_cmd_t                   cmd_le_gap_set_adv_data;
    struct gecko_msg_le_gap_set_adv_data_rsp_t                   rsp_le_gap_set_adv_data;
    struct gecko_msg_le_gap_set_adv_timeout_cmd_t                cmd_le_gap_set_adv_timeout;
    struct gecko_msg_le_gap_set_adv_timeout_rsp_t                rsp_le_gap_set_adv_timeout;
    struct gecko_msg_le_gap_set_conn_phy_cmd_t                   cmd_le_gap_set_conn_phy;
    struct gecko_msg_le_gap_set_conn_phy_rsp_t                   rsp_le_gap_set_conn_phy;
    struct gecko_msg_le_gap_bt5_set_mode_cmd_t                   cmd_le_gap_bt5_set_mode;
    struct gecko_msg_le_gap_bt5_set_mode_rsp_t                   rsp_le_gap_bt5_set_mode;
    struct gecko_msg_le_gap_bt5_set_adv_parameters_cmd_t         cmd_le_gap_bt5_set_adv_parameters;
    struct gecko_msg_le_gap_bt5_set_adv_parameters_rsp_t         rsp_le_gap_bt5_set_adv_parameters;
    struct gecko_msg_le_gap_bt5_set_adv_data_cmd_t               cmd_le_gap_bt5_set_adv_data;
    struct gecko_msg_le_gap_bt5_set_adv_data_rsp_t               rsp_le_gap_bt5_set_adv_data;
    struct gecko_msg_le_gap_set_privacy_mode_cmd_t               cmd_le_gap_set_privacy_mode;
    struct gecko_msg_le_gap_set_privacy_mode_rsp_t               rsp_le_gap_set_privacy_mode;
    struct gecko_msg_le_gap_set_advertise_timing_cmd_t           cmd_le_gap_set_advertise_timing;
    struct gecko_msg_le_gap_set_advertise_timing_rsp_t           rsp_le_gap_set_advertise_timing;
    struct gecko_msg_le_gap_set_advertise_channel_map_cmd_t      cmd_le_gap_set_advertise_channel_map;
    struct gecko_msg_le_gap_set_advertise_channel_map_rsp_t      rsp_le_gap_set_advertise_channel_map;
    struct gecko_msg_le_gap_set_advertise_report_scan_request_cmd_t cmd_le_gap_set_advertise_report_scan_request;
    struct gecko_msg_le_gap_set_advertise_report_scan_request_rsp_t rsp_le_gap_set_advertise_report_scan_request;
    struct gecko_msg_le_gap_set_advertise_phy_cmd_t              cmd_le_gap_set_advertise_phy;
    struct gecko_msg_le_gap_set_advertise_phy_rsp_t              rsp_le_gap_set_advertise_phy;
    struct gecko_msg_le_gap_set_advertise_configuration_cmd_t    cmd_le_gap_set_advertise_configuration;
    struct gecko_msg_le_gap_set_advertise_configuration_rsp_t    rsp_le_gap_set_advertise_configuration;
    struct gecko_msg_le_gap_clear_advertise_configuration_cmd_t  cmd_le_gap_clear_advertise_configuration;
    struct gecko_msg_le_gap_clear_advertise_configuration_rsp_t  rsp_le_gap_clear_advertise_configuration;
    struct gecko_msg_le_gap_start_advertising_cmd_t              cmd_le_gap_start_advertising;
    struct gecko_msg_le_gap_start_advertising_rsp_t              rsp_le_gap_start_advertising;
    struct gecko_msg_le_gap_stop_advertising_cmd_t               cmd_le_gap_stop_advertising;
    struct gecko_msg_le_gap_stop_advertising_rsp_t               rsp_le_gap_stop_advertising;
    struct gecko_msg_le_gap_set_discovery_timing_cmd_t           cmd_le_gap_set_discovery_timing;
    struct gecko_msg_le_gap_set_discovery_timing_rsp_t           rsp_le_gap_set_discovery_timing;
    struct gecko_msg_le_gap_set_discovery_type_cmd_t             cmd_le_gap_set_discovery_type;
    struct gecko_msg_le_gap_set_discovery_type_rsp_t             rsp_le_gap_set_discovery_type;
    struct gecko_msg_le_gap_start_discovery_cmd_t                cmd_le_gap_start_discovery;
    struct gecko_msg_le_gap_start_discovery_rsp_t                rsp_le_gap_start_discovery;
    struct gecko_msg_le_gap_set_data_channel_classification_cmd_t cmd_le_gap_set_data_channel_classification;
    struct gecko_msg_le_gap_set_data_channel_classification_rsp_t rsp_le_gap_set_data_channel_classification;
    struct gecko_msg_le_gap_connect_cmd_t                        cmd_le_gap_connect;
    struct gecko_msg_le_gap_connect_rsp_t                        rsp_le_gap_connect;
    struct gecko_msg_le_gap_set_advertise_tx_power_cmd_t         cmd_le_gap_set_advertise_tx_power;
    struct gecko_msg_le_gap_set_advertise_tx_power_rsp_t         rsp_le_gap_set_advertise_tx_power;
    struct gecko_msg_le_gap_set_discovery_extended_scan_response_cmd_t cmd_le_gap_set_discovery_extended_scan_response;
    struct gecko_msg_le_gap_set_discovery_extended_scan_response_rsp_t rsp_le_gap_set_discovery_extended_scan_response;
    struct gecko_msg_le_gap_start_periodic_advertising_cmd_t     cmd_le_gap_start_periodic_advertising;
    struct gecko_msg_le_gap_start_periodic_advertising_rsp_t     rsp_le_gap_start_periodic_advertising;
    struct gecko_msg_le_gap_stop_periodic_advertising_cmd_t      cmd_le_gap_stop_periodic_advertising;
    struct gecko_msg_le_gap_stop_periodic_advertising_rsp_t      rsp_le_gap_stop_periodic_advertising;
    struct gecko_msg_le_gap_set_long_advertising_data_cmd_t      cmd_le_gap_set_long_advertising_data;
    struct gecko_msg_le_gap_set_long_advertising_data_rsp_t      rsp_le_gap_set_long_advertising_data;
    struct gecko_msg_le_gap_enable_whitelisting_cmd_t            cmd_le_gap_enable_whitelisting;
    struct gecko_msg_le_gap_enable_whitelisting_rsp_t            rsp_le_gap_enable_whitelisting;
    struct gecko_msg_le_gap_set_conn_timing_parameters_cmd_t     cmd_le_gap_set_conn_timing_parameters;
    struct gecko_msg_le_gap_set_conn_timing_parameters_rsp_t     rsp_le_gap_set_conn_timing_parameters;
    struct gecko_msg_le_gap_set_advertise_random_address_cmd_t   cmd_le_gap_set_advertise_random_address;
    struct gecko_msg_le_gap_set_advertise_random_address_rsp_t   rsp_le_gap_set_advertise_random_address;
    struct gecko_msg_le_gap_clear_advertise_random_address_cmd_t cmd_le_gap_clear_advertise_random_address;
    struct gecko_msg_le_gap_clear_advertise_random_address_rsp_t rsp_le_gap_clear_advertise_random_address;
    struct gecko_msg_le_gap_scan_response_evt_t                  evt_le_gap_scan_response;
    struct gecko_msg_le_gap_adv_timeout_evt_t                    evt_le_gap_adv_timeout;
    struct gecko_msg_le_gap_scan_request_evt_t                   evt_le_gap_scan_request;
    struct gecko_msg_le_gap_extended_scan_response_evt_t         evt_le_gap_extended_scan_response;
    struct gecko_msg_le_gap_periodic_advertising_status_evt_t    evt_le_gap_periodic_advertising_status;
    struct gecko_msg_sync_open_cmd_t                             cmd_sync_open;
    struct gecko_msg_sync_open_rsp_t                             rsp_sync_open;
    struct gecko_msg_sync_close_cmd_t                            cmd_sync_close;
    struct gecko_msg_sync_close_rsp_t                            rsp_sync_close;
    struct gecko_msg_sync_opened_evt_t                           evt_sync_opened;
    struct gecko_msg_sync_closed_evt_t                           evt_sync_closed;
    struct gecko_msg_sync_data_evt_t                             evt_sync_data;
    struct gecko_msg_le_connection_set_parameters_cmd_t          cmd_le_connection_set_parameters;
    struct gecko_msg_le_connection_set_parameters_rsp_t          rsp_le_connection_set_parameters;
    struct gecko_msg_le_connection_get_rssi_cmd_t                cmd_le_connection_get_rssi;
    struct gecko_msg_le_connection_get_rssi_rsp_t                rsp_le_connection_get_rssi;
    struct gecko_msg_le_connection_disable_slave_latency_cmd_t   cmd_le_connection_disable_slave_latency;
    struct gecko_msg_le_connection_disable_slave_latency_rsp_t   rsp_le_connection_disable_slave_latency;
    struct gecko_msg_le_connection_set_phy_cmd_t                 cmd_le_connection_set_phy;
    struct gecko_msg_le_connection_set_phy_rsp_t                 rsp_le_connection_set_phy;
    struct gecko_msg_le_connection_close_cmd_t                   cmd_le_connection_close;
    struct gecko_msg_le_connection_close_rsp_t                   rsp_le_connection_close;
    struct gecko_msg_le_connection_set_timing_parameters_cmd_t   cmd_le_connection_set_timing_parameters;
    struct gecko_msg_le_connection_set_timing_parameters_rsp_t   rsp_le_connection_set_timing_parameters;
    struct gecko_msg_le_connection_read_channel_map_cmd_t        cmd_le_connection_read_channel_map;
    struct gecko_msg_le_connection_read_channel_map_rsp_t        rsp_le_connection_read_channel_map;
    struct gecko_msg_le_connection_set_preferred_phy_cmd_t       cmd_le_connection_set_preferred_phy;
    struct gecko_msg_le_connection_set_preferred_phy_rsp_t       rsp_le_connection_set_preferred_phy;
    struct gecko_msg_le_connection_opened_evt_t                  evt_le_connection_opened;
    struct gecko_msg_le_connection_closed_evt_t                  evt_le_connection_closed;
    struct gecko_msg_le_connection_parameters_evt_t              evt_le_connection_parameters;
    struct gecko_msg_le_connection_rssi_evt_t                    evt_le_connection_rssi;
    struct gecko_msg_le_connection_phy_status_evt_t              evt_le_connection_phy_status;
    struct gecko_msg_gatt_set_max_mtu_cmd_t                      cmd_gatt_set_max_mtu;
    struct gecko_msg_gatt_set_max_mtu_rsp_t                      rsp_gatt_set_max_mtu;
    struct gecko_msg_gatt_discover_primary_services_cmd_t        cmd_gatt_discover_primary_services;
    struct gecko_msg_gatt_discover_primary_services_rsp_t        rsp_gatt_discover_primary_services;
    struct gecko_msg_gatt_discover_primary_services_by_uuid_cmd_t cmd_gatt_discover_primary_services_by_uuid;
    struct gecko_msg_gatt_discover_primary_services_by_uuid_rsp_t rsp_gatt_discover_primary_services_by_uuid;
    struct gecko_msg_gatt_discover_characteristics_cmd_t         cmd_gatt_discover_characteristics;
    struct gecko_msg_gatt_discover_characteristics_rsp_t         rsp_gatt_discover_characteristics;
    struct gecko_msg_gatt_discover_characteristics_by_uuid_cmd_t cmd_gatt_discover_characteristics_by_uuid;
    struct gecko_msg_gatt_discover_characteristics_by_uuid_rsp_t rsp_gatt_discover_characteristics_by_uuid;
    struct gecko_msg_gatt_set_characteristic_notification_cmd_t  cmd_gatt_set_characteristic_notification;
    struct gecko_msg_gatt_set_characteristic_notification_rsp_t  rsp_gatt_set_characteristic_notification;
    struct gecko_msg_gatt_discover_descriptors_cmd_t             cmd_gatt_discover_descriptors;
    struct gecko_msg_gatt_discover_descriptors_rsp_t             rsp_gatt_discover_descriptors;
    struct gecko_msg_gatt_read_characteristic_value_cmd_t        cmd_gatt_read_characteristic_value;
    struct gecko_msg_gatt_read_characteristic_value_rsp_t        rsp_gatt_read_characteristic_value;
    struct gecko_msg_gatt_read_characteristic_value_by_uuid_cmd_t cmd_gatt_read_characteristic_value_by_uuid;
    struct gecko_msg_gatt_read_characteristic_value_by_uuid_rsp_t rsp_gatt_read_characteristic_value_by_uuid;
    struct gecko_msg_gatt_write_characteristic_value_cmd_t       cmd_gatt_write_characteristic_value;
    struct gecko_msg_gatt_write_characteristic_value_rsp_t       rsp_gatt_write_characteristic_value;
    struct gecko_msg_gatt_write_characteristic_value_without_response_cmd_t cmd_gatt_write_characteristic_value_without_response;
    struct gecko_msg_gatt_write_characteristic_value_without_response_rsp_t rsp_gatt_write_characteristic_value_without_response;
    struct gecko_msg_gatt_prepare_characteristic_value_write_cmd_t cmd_gatt_prepare_characteristic_value_write;
    struct gecko_msg_gatt_prepare_characteristic_value_write_rsp_t rsp_gatt_prepare_characteristic_value_write;
    struct gecko_msg_gatt_execute_characteristic_value_write_cmd_t cmd_gatt_execute_characteristic_value_write;
    struct gecko_msg_gatt_execute_characteristic_value_write_rsp_t rsp_gatt_execute_characteristic_value_write;
    struct gecko_msg_gatt_send_characteristic_confirmation_cmd_t cmd_gatt_send_characteristic_confirmation;
    struct gecko_msg_gatt_send_characteristic_confirmation_rsp_t rsp_gatt_send_characteristic_confirmation;
    struct gecko_msg_gatt_read_descriptor_value_cmd_t            cmd_gatt_read_descriptor_value;
    struct gecko_msg_gatt_read_descriptor_value_rsp_t            rsp_gatt_read_descriptor_value;
    struct gecko_msg_gatt_write_descriptor_value_cmd_t           cmd_gatt_write_descriptor_value;
    struct gecko_msg_gatt_write_descriptor_value_rsp_t           rsp_gatt_write_descriptor_value;
    struct gecko_msg_gatt_find_included_services_cmd_t           cmd_gatt_find_included_services;
    struct gecko_msg_gatt_find_included_services_rsp_t           rsp_gatt_find_included_services;
    struct gecko_msg_gatt_read_multiple_characteristic_values_cmd_t cmd_gatt_read_multiple_characteristic_values;
    struct gecko_msg_gatt_read_multiple_characteristic_values_rsp_t rsp_gatt_read_multiple_characteristic_values;
    struct gecko_msg_gatt_read_characteristic_value_from_offset_cmd_t cmd_gatt_read_characteristic_value_from_offset;
    struct gecko_msg_gatt_read_characteristic_value_from_offset_rsp_t rsp_gatt_read_characteristic_value_from_offset;
    struct gecko_msg_gatt_prepare_characteristic_value_reliable_write_cmd_t cmd_gatt_prepare_characteristic_value_reliable_write;
    struct gecko_msg_gatt_prepare_characteristic_value_reliable_write_rsp_t rsp_gatt_prepare_characteristic_value_reliable_write;
    struct gecko_msg_gatt_mtu_exchanged_evt_t                    evt_gatt_mtu_exchanged;
    struct gecko_msg_gatt_service_evt_t                          evt_gatt_service;
    struct gecko_msg_gatt_characteristic_evt_t                   evt_gatt_characteristic;
    struct gecko_msg_gatt_descriptor_evt_t                       evt_gatt_descriptor;
    struct gecko_msg_gatt_characteristic_value_evt_t             evt_gatt_characteristic_value;
    struct gecko_msg_gatt_descriptor_value_evt_t                 evt_gatt_descriptor_value;
    struct gecko_msg_gatt_procedure_completed_evt_t              evt_gatt_procedure_completed;
    struct gecko_msg_gatt_server_read_attribute_value_cmd_t      cmd_gatt_server_read_attribute_value;
    struct gecko_msg_gatt_server_read_attribute_value_rsp_t      rsp_gatt_server_read_attribute_value;
    struct gecko_msg_gatt_server_read_attribute_type_cmd_t       cmd_gatt_server_read_attribute_type;
    struct gecko_msg_gatt_server_read_attribute_type_rsp_t       rsp_gatt_server_read_attribute_type;
    struct gecko_msg_gatt_server_write_attribute_value_cmd_t     cmd_gatt_server_write_attribute_value;
    struct gecko_msg_gatt_server_write_attribute_value_rsp_t     rsp_gatt_server_write_attribute_value;
    struct gecko_msg_gatt_server_send_user_read_response_cmd_t   cmd_gatt_server_send_user_read_response;
    struct gecko_msg_gatt_server_send_user_read_response_rsp_t   rsp_gatt_server_send_user_read_response;
    struct gecko_msg_gatt_server_send_user_write_response_cmd_t  cmd_gatt_server_send_user_write_response;
    struct gecko_msg_gatt_server_send_user_write_response_rsp_t  rsp_gatt_server_send_user_write_response;
    struct gecko_msg_gatt_server_send_characteristic_notification_cmd_t cmd_gatt_server_send_characteristic_notification;
    struct gecko_msg_gatt_server_send_characteristic_notification_rsp_t rsp_gatt_server_send_characteristic_notification;
    struct gecko_msg_gatt_server_find_attribute_cmd_t            cmd_gatt_server_find_attribute;
    struct gecko_msg_gatt_server_find_attribute_rsp_t            rsp_gatt_server_find_attribute;
    struct gecko_msg_gatt_server_set_capabilities_cmd_t          cmd_gatt_server_set_capabilities;
    struct gecko_msg_gatt_server_set_capabilities_rsp_t          rsp_gatt_server_set_capabilities;
    struct gecko_msg_gatt_server_set_max_mtu_cmd_t               cmd_gatt_server_set_max_mtu;
    struct gecko_msg_gatt_server_set_max_mtu_rsp_t               rsp_gatt_server_set_max_mtu;
    struct gecko_msg_gatt_server_get_mtu_cmd_t                   cmd_gatt_server_get_mtu;
    struct gecko_msg_gatt_server_get_mtu_rsp_t                   rsp_gatt_server_get_mtu;
    struct gecko_msg_gatt_server_enable_capabilities_cmd_t       cmd_gatt_server_enable_capabilities;
    struct gecko_msg_gatt_server_enable_capabilities_rsp_t       rsp_gatt_server_enable_capabilities;
    struct gecko_msg_gatt_server_disable_capabilities_cmd_t      cmd_gatt_server_disable_capabilities;
    struct gecko_msg_gatt_server_disable_capabilities_rsp_t      rsp_gatt_server_disable_capabilities;
    struct gecko_msg_gatt_server_get_enabled_capabilities_rsp_t  rsp_gatt_server_get_enabled_capabilities;
    struct gecko_msg_gatt_server_attribute_value_evt_t           evt_gatt_server_attribute_value;
    struct gecko_msg_gatt_server_user_read_request_evt_t         evt_gatt_server_user_read_request;
    struct gecko_msg_gatt_server_user_write_request_evt_t        evt_gatt_server_user_write_request;
    struct gecko_msg_gatt_server_characteristic_status_evt_t     evt_gatt_server_characteristic_status;
    struct gecko_msg_gatt_server_execute_write_completed_evt_t   evt_gatt_server_execute_write_completed;
    struct gecko_msg_hardware_set_soft_timer_cmd_t               cmd_hardware_set_soft_timer;
    struct gecko_msg_hardware_set_soft_timer_rsp_t               rsp_hardware_set_soft_timer;
    struct gecko_msg_hardware_get_time_rsp_t                     rsp_hardware_get_time;
    struct gecko_msg_hardware_set_lazy_soft_timer_cmd_t          cmd_hardware_set_lazy_soft_timer;
    struct gecko_msg_hardware_set_lazy_soft_timer_rsp_t          rsp_hardware_set_lazy_soft_timer;
    struct gecko_msg_hardware_soft_timer_evt_t                   evt_hardware_soft_timer;
    struct gecko_msg_flash_ps_erase_all_rsp_t                    rsp_flash_ps_erase_all;
    struct gecko_msg_flash_ps_save_cmd_t                         cmd_flash_ps_save;
    struct gecko_msg_flash_ps_save_rsp_t                         rsp_flash_ps_save;
    struct gecko_msg_flash_ps_load_cmd_t                         cmd_flash_ps_load;
    struct gecko_msg_flash_ps_load_rsp_t                         rsp_flash_ps_load;
    struct gecko_msg_flash_ps_erase_cmd_t                        cmd_flash_ps_erase;
    struct gecko_msg_flash_ps_erase_rsp_t                        rsp_flash_ps_erase;
    struct gecko_msg_test_dtm_tx_cmd_t                           cmd_test_dtm_tx;
    struct gecko_msg_test_dtm_tx_rsp_t                           rsp_test_dtm_tx;
    struct gecko_msg_test_dtm_rx_cmd_t                           cmd_test_dtm_rx;
    struct gecko_msg_test_dtm_rx_rsp_t                           rsp_test_dtm_rx;
    struct gecko_msg_test_dtm_end_rsp_t                          rsp_test_dtm_end;
    struct gecko_msg_test_dtm_completed_evt_t                    evt_test_dtm_completed;
    struct gecko_msg_sm_set_bondable_mode_cmd_t                  cmd_sm_set_bondable_mode;
    struct gecko_msg_sm_set_bondable_mode_rsp_t                  rsp_sm_set_bondable_mode;
    struct gecko_msg_sm_configure_cmd_t                          cmd_sm_configure;
    struct gecko_msg_sm_configure_rsp_t                          rsp_sm_configure;
    struct gecko_msg_sm_store_bonding_configuration_cmd_t        cmd_sm_store_bonding_configuration;
    struct gecko_msg_sm_store_bonding_configuration_rsp_t        rsp_sm_store_bonding_configuration;
    struct gecko_msg_sm_increase_security_cmd_t                  cmd_sm_increase_security;
    struct gecko_msg_sm_increase_security_rsp_t                  rsp_sm_increase_security;
    struct gecko_msg_sm_delete_bonding_cmd_t                     cmd_sm_delete_bonding;
    struct gecko_msg_sm_delete_bonding_rsp_t                     rsp_sm_delete_bonding;
    struct gecko_msg_sm_delete_bondings_rsp_t                    rsp_sm_delete_bondings;
    struct gecko_msg_sm_enter_passkey_cmd_t                      cmd_sm_enter_passkey;
    struct gecko_msg_sm_enter_passkey_rsp_t                      rsp_sm_enter_passkey;
    struct gecko_msg_sm_passkey_confirm_cmd_t                    cmd_sm_passkey_confirm;
    struct gecko_msg_sm_passkey_confirm_rsp_t                    rsp_sm_passkey_confirm;
    struct gecko_msg_sm_set_oob_data_cmd_t                       cmd_sm_set_oob_data;
    struct gecko_msg_sm_set_oob_data_rsp_t                       rsp_sm_set_oob_data;
    struct gecko_msg_sm_list_all_bondings_rsp_t                  rsp_sm_list_all_bondings;
    struct gecko_msg_sm_bonding_confirm_cmd_t                    cmd_sm_bonding_confirm;
    struct gecko_msg_sm_bonding_confirm_rsp_t                    rsp_sm_bonding_confirm;
    struct gecko_msg_sm_set_debug_mode_rsp_t                     rsp_sm_set_debug_mode;
    struct gecko_msg_sm_set_passkey_cmd_t                        cmd_sm_set_passkey;
    struct gecko_msg_sm_set_passkey_rsp_t                        rsp_sm_set_passkey;
    struct gecko_msg_sm_use_sc_oob_cmd_t                         cmd_sm_use_sc_oob;
    struct gecko_msg_sm_use_sc_oob_rsp_t                         rsp_sm_use_sc_oob;
    struct gecko_msg_sm_set_sc_remote_oob_data_cmd_t             cmd_sm_set_sc_remote_oob_data;
    struct gecko_msg_sm_set_sc_remote_oob_data_rsp_t             rsp_sm_set_sc_remote_oob_data;
    struct gecko_msg_sm_add_to_whitelist_cmd_t                   cmd_sm_add_to_whitelist;
    struct gecko_msg_sm_add_to_whitelist_rsp_t                   rsp_sm_add_to_whitelist;
    struct gecko_msg_sm_set_minimum_key_size_cmd_t               cmd_sm_set_minimum_key_size;
    struct gecko_msg_sm_set_minimum_key_size_rsp_t               rsp_sm_set_minimum_key_size;
    struct gecko_msg_sm_passkey_display_evt_t                    evt_sm_passkey_display;
    struct gecko_msg_sm_passkey_request_evt_t                    evt_sm_passkey_request;
    struct gecko_msg_sm_confirm_passkey_evt_t                    evt_sm_confirm_passkey;
    struct gecko_msg_sm_bonded_evt_t                             evt_sm_bonded;
    struct gecko_msg_sm_bonding_failed_evt_t                     evt_sm_bonding_failed;
    struct gecko_msg_sm_list_bonding_entry_evt_t                 evt_sm_list_bonding_entry;
    struct gecko_msg_sm_confirm_bonding_evt_t                    evt_sm_confirm_bonding;
    struct gecko_msg_homekit_configure_cmd_t                     cmd_homekit_configure;
    struct gecko_msg_homekit_configure_rsp_t                     rsp_homekit_configure;
    struct gecko_msg_homekit_advertise_cmd_t                     cmd_homekit_advertise;
    struct gecko_msg_homekit_advertise_rsp_t                     rsp_homekit_advertise;
    struct gecko_msg_homekit_delete_pairings_rsp_t               rsp_homekit_delete_pairings;
    struct gecko_msg_homekit_check_authcp_rsp_t                  rsp_homekit_check_authcp;
    struct gecko_msg_homekit_get_pairing_id_cmd_t                cmd_homekit_get_pairing_id;
    struct gecko_msg_homekit_get_pairing_id_rsp_t                rsp_homekit_get_pairing_id;
    struct gecko_msg_homekit_send_write_response_cmd_t           cmd_homekit_send_write_response;
    struct gecko_msg_homekit_send_write_response_rsp_t           rsp_homekit_send_write_response;
    struct gecko_msg_homekit_send_read_response_cmd_t            cmd_homekit_send_read_response;
    struct gecko_msg_homekit_send_read_response_rsp_t            rsp_homekit_send_read_response;
    struct gecko_msg_homekit_gsn_action_cmd_t                    cmd_homekit_gsn_action;
    struct gecko_msg_homekit_gsn_action_rsp_t                    rsp_homekit_gsn_action;
    struct gecko_msg_homekit_event_notification_cmd_t            cmd_homekit_event_notification;
    struct gecko_msg_homekit_event_notification_rsp_t            rsp_homekit_event_notification;
    struct gecko_msg_homekit_broadcast_action_cmd_t              cmd_homekit_broadcast_action;
    struct gecko_msg_homekit_broadcast_action_rsp_t              rsp_homekit_broadcast_action;
    struct gecko_msg_homekit_configure_product_data_cmd_t        cmd_homekit_configure_product_data;
    struct gecko_msg_homekit_configure_product_data_rsp_t        rsp_homekit_configure_product_data;
    struct gecko_msg_homekit_setupcode_display_evt_t             evt_homekit_setupcode_display;
    struct gecko_msg_homekit_paired_evt_t                        evt_homekit_paired;
    struct gecko_msg_homekit_pair_verified_evt_t                 evt_homekit_pair_verified;
    struct gecko_msg_homekit_connection_opened_evt_t             evt_homekit_connection_opened;
    struct gecko_msg_homekit_connection_closed_evt_t             evt_homekit_connection_closed;
    struct gecko_msg_homekit_identify_evt_t                      evt_homekit_identify;
    struct gecko_msg_homekit_write_request_evt_t                 evt_homekit_write_request;
    struct gecko_msg_homekit_read_request_evt_t                  evt_homekit_read_request;
    struct gecko_msg_homekit_disconnection_required_evt_t        evt_homekit_disconnection_required;
    struct gecko_msg_homekit_pairing_removed_evt_t               evt_homekit_pairing_removed;
    struct gecko_msg_homekit_setuppayload_display_evt_t          evt_homekit_setuppayload_display;
    struct gecko_msg_coex_set_options_cmd_t                      cmd_coex_set_options;
    struct gecko_msg_coex_set_options_rsp_t                      rsp_coex_set_options;
    struct gecko_msg_coex_get_counters_cmd_t                     cmd_coex_get_counters;
    struct gecko_msg_coex_get_counters_rsp_t                     rsp_coex_get_counters;
    struct gecko_msg_coex_set_parameters_cmd_t                   cmd_coex_set_parameters;
    struct gecko_msg_coex_set_parameters_rsp_t                   rsp_coex_set_parameters;
    struct gecko_msg_coex_set_directional_priority_pulse_cmd_t   cmd_coex_set_directional_priority_pulse;
    struct gecko_msg_coex_set_directional_priority_pulse_rsp_t   rsp_coex_set_directional_priority_pulse;
    struct gecko_msg_l2cap_coc_send_connection_request_cmd_t     cmd_l2cap_coc_send_connection_request;
    struct gecko_msg_l2cap_coc_send_connection_request_rsp_t     rsp_l2cap_coc_send_connection_request;
    struct gecko_msg_l2cap_coc_send_connection_response_cmd_t    cmd_l2cap_coc_send_connection_response;
    struct gecko_msg_l2cap_coc_send_connection_response_rsp_t    rsp_l2cap_coc_send_connection_response;
    struct gecko_msg_l2cap_coc_send_le_flow_control_credit_cmd_t cmd_l2cap_coc_send_le_flow_control_credit;
    struct gecko_msg_l2cap_coc_send_le_flow_control_credit_rsp_t rsp_l2cap_coc_send_le_flow_control_credit;
    struct gecko_msg_l2cap_coc_send_disconnection_request_cmd_t  cmd_l2cap_coc_send_disconnection_request;
    struct gecko_msg_l2cap_coc_send_disconnection_request_rsp_t  rsp_l2cap_coc_send_disconnection_request;
    struct gecko_msg_l2cap_coc_send_data_cmd_t                   cmd_l2cap_coc_send_data;
    struct gecko_msg_l2cap_coc_send_data_rsp_t                   rsp_l2cap_coc_send_data;
    struct gecko_msg_l2cap_coc_connection_request_evt_t          evt_l2cap_coc_connection_request;
    struct gecko_msg_l2cap_coc_connection_response_evt_t         evt_l2cap_coc_connection_response;
    struct gecko_msg_l2cap_coc_le_flow_control_credit_evt_t      evt_l2cap_coc_le_flow_control_credit;
    struct gecko_msg_l2cap_coc_channel_disconnected_evt_t        evt_l2cap_coc_channel_disconnected;
    struct gecko_msg_l2cap_coc_data_evt_t                        evt_l2cap_coc_data;
    struct gecko_msg_l2cap_command_rejected_evt_t                evt_l2cap_command_rejected;
    struct gecko_msg_cte_transmitter_enable_connection_cte_cmd_t cmd_cte_transmitter_enable_connection_cte;
    struct gecko_msg_cte_transmitter_enable_connection_cte_rsp_t rsp_cte_transmitter_enable_connection_cte;
    struct gecko_msg_cte_transmitter_disable_connection_cte_cmd_t cmd_cte_transmitter_disable_connection_cte;
    struct gecko_msg_cte_transmitter_disable_connection_cte_rsp_t rsp_cte_transmitter_disable_connection_cte;
    struct gecko_msg_cte_transmitter_enable_connectionless_cte_cmd_t cmd_cte_transmitter_enable_connectionless_cte;
    struct gecko_msg_cte_transmitter_enable_connectionless_cte_rsp_t rsp_cte_transmitter_enable_connectionless_cte;
    struct gecko_msg_cte_transmitter_disable_connectionless_cte_cmd_t cmd_cte_transmitter_disable_connectionless_cte;
    struct gecko_msg_cte_transmitter_disable_connectionless_cte_rsp_t rsp_cte_transmitter_disable_connectionless_cte;
    struct gecko_msg_cte_transmitter_set_dtm_parameters_cmd_t    cmd_cte_transmitter_set_dtm_parameters;
    struct gecko_msg_cte_transmitter_set_dtm_parameters_rsp_t    rsp_cte_transmitter_set_dtm_parameters;
    struct gecko_msg_cte_transmitter_clear_dtm_parameters_rsp_t  rsp_cte_transmitter_clear_dtm_parameters;
    struct gecko_msg_cte_transmitter_enable_silabs_cte_cmd_t     cmd_cte_transmitter_enable_silabs_cte;
    struct gecko_msg_cte_transmitter_enable_silabs_cte_rsp_t     rsp_cte_transmitter_enable_silabs_cte;
    struct gecko_msg_cte_transmitter_disable_silabs_cte_cmd_t    cmd_cte_transmitter_disable_silabs_cte;
    struct gecko_msg_cte_transmitter_disable_silabs_cte_rsp_t    rsp_cte_transmitter_disable_silabs_cte;
    struct gecko_msg_cte_receiver_configure_cmd_t                cmd_cte_receiver_configure;
    struct gecko_msg_cte_receiver_configure_rsp_t                rsp_cte_receiver_configure;
    struct gecko_msg_cte_receiver_enable_connection_cte_cmd_t    cmd_cte_receiver_enable_connection_cte;
    struct gecko_msg_cte_receiver_enable_connection_cte_rsp_t    rsp_cte_receiver_enable_connection_cte;
    struct gecko_msg_cte_receiver_disable_connection_cte_cmd_t   cmd_cte_receiver_disable_connection_cte;
    struct gecko_msg_cte_receiver_disable_connection_cte_rsp_t   rsp_cte_receiver_disable_connection_cte;
    struct gecko_msg_cte_receiver_enable_connectionless_cte_cmd_t cmd_cte_receiver_enable_connectionless_cte;
    struct gecko_msg_cte_receiver_enable_connectionless_cte_rsp_t rsp_cte_receiver_enable_connectionless_cte;
    struct gecko_msg_cte_receiver_disable_connectionless_cte_cmd_t cmd_cte_receiver_disable_connectionless_cte;
    struct gecko_msg_cte_receiver_disable_connectionless_cte_rsp_t rsp_cte_receiver_disable_connectionless_cte;
    struct gecko_msg_cte_receiver_set_dtm_parameters_cmd_t       cmd_cte_receiver_set_dtm_parameters;
    struct gecko_msg_cte_receiver_set_dtm_parameters_rsp_t       rsp_cte_receiver_set_dtm_parameters;
    struct gecko_msg_cte_receiver_clear_dtm_parameters_rsp_t     rsp_cte_receiver_clear_dtm_parameters;
    struct gecko_msg_cte_receiver_enable_silabs_cte_cmd_t        cmd_cte_receiver_enable_silabs_cte;
    struct gecko_msg_cte_receiver_enable_silabs_cte_rsp_t        rsp_cte_receiver_enable_silabs_cte;
    struct gecko_msg_cte_receiver_disable_silabs_cte_rsp_t       rsp_cte_receiver_disable_silabs_cte;
    struct gecko_msg_cte_receiver_connection_iq_report_evt_t     evt_cte_receiver_connection_iq_report;
    struct gecko_msg_cte_receiver_connectionless_iq_report_evt_t evt_cte_receiver_connectionless_iq_report;
    struct gecko_msg_cte_receiver_dtm_iq_report_evt_t            evt_cte_receiver_dtm_iq_report;
    struct gecko_msg_cte_receiver_silabs_iq_report_evt_t         evt_cte_receiver_silabs_iq_report;
    struct gecko_msg_user_message_to_target_cmd_t                cmd_user_message_to_target;
    struct gecko_msg_user_message_to_target_rsp_t                rsp_user_message_to_target;
    struct gecko_msg_user_message_to_host_evt_t                  evt_user_message_to_host;

    uint8 payload[BGLIB_MSG_MAX_PAYLOAD];
}data;

});


void gecko_handle_command(uint32_t,void*);
void gecko_handle_command_noresponse(uint32_t,void*);

extern struct gecko_cmd_packet*  gecko_cmd_msg;
extern struct gecko_cmd_packet*  gecko_rsp_msg;

/** 
*
* gecko_cmd_dfu_reset
*
* Reset the system. The command does not have a response but it triggers one of
* the boot events (normal reset or boot to DFU mode) after re-boot. 
*
* @param dfu   Boot mode:
*  
*      0: Normal reset  
*      1: Boot to UART DFU mode  
*      2: Boot to OTA DFU mode
*
* Events generated
*
* gecko_evt_system_boot - Sent after the device has booted in normal mode
* gecko_evt_dfu_boot - Sent after the device has booted in UART DFU mode
*
**/

static inline void* gecko_cmd_dfu_reset(uint8 dfu)
{
    
    gecko_cmd_msg->data.cmd_dfu_reset.dfu=dfu;
    gecko_cmd_msg->header=(gecko_cmd_dfu_reset_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command_noresponse(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    return 0;
}

/** 
*
* gecko_cmd_dfu_flash_set_address
*
* After re-booting the local device in DFU mode, this command defines the
* starting address on the flash where the new firmware will be written. 
*
* @param address   The offset in the flash where the new firmware is uploaded to. Always use the
*  value 0x00000000.
*
**/

static inline struct gecko_msg_dfu_flash_set_address_rsp_t* gecko_cmd_dfu_flash_set_address(uint32 address)
{
    
    gecko_cmd_msg->data.cmd_dfu_flash_set_address.address=address;
    gecko_cmd_msg->header=(gecko_cmd_dfu_flash_set_address_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_dfu_flash_set_address;
}

/** 
*
* gecko_cmd_dfu_flash_upload
*
* Upload the whole firmware image file into the Bluetooth device. The passed
* data length must be a multiple of 4 bytes. Because the BGAPI command payload
* size is limited, multiple commands need to be issued one after the other until
* the whole .bin firmware image file is uploaded to the device. After each
* command, the next address of the flash sector in memory to write to is
* automatically updated by the bootloader. 
*
* @param data_len   Array length
* @param data_data   An array of data which will be written onto the flash.
*
**/

static inline struct gecko_msg_dfu_flash_upload_rsp_t* gecko_cmd_dfu_flash_upload(uint8 data_len, const uint8* data_data)
{
    if ((uint16_t)data_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_dfu_flash_upload.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_dfu_flash_upload;
    }

    
    gecko_cmd_msg->data.cmd_dfu_flash_upload.data.len=data_len;
    memcpy(gecko_cmd_msg->data.cmd_dfu_flash_upload.data.data,data_data,data_len);
    gecko_cmd_msg->header=(gecko_cmd_dfu_flash_upload_id+(((1+data_len)&0xff)<<8)+(((1+data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_dfu_flash_upload;
}

/** 
*
* gecko_cmd_dfu_flash_upload_finish
*
* Inform the device that the DFU file is fully uploaded. To return the device
* back to normal mode, issue the command DFU Reset . 
*
*
**/

static inline struct gecko_msg_dfu_flash_upload_finish_rsp_t* gecko_cmd_dfu_flash_upload_finish()
{
    
    gecko_cmd_msg->header=(gecko_cmd_dfu_flash_upload_finish_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_dfu_flash_upload_finish;
}

/** 
*
* gecko_cmd_system_hello
*
* Verify whether the communication between the host and the device is
* functional. 
*
*
**/

static inline struct gecko_msg_system_hello_rsp_t* gecko_cmd_system_hello()
{
    
    gecko_cmd_msg->header=(gecko_cmd_system_hello_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_hello;
}

/** 
*
* gecko_cmd_system_reset
*
* Reset the system. The command does not have a response but it triggers one of
* the boot events (normal reset or boot to DFU mode) depending on the selected
* boot mode. 
*
* @param dfu   Boot mode:
*  
*      0: Normal reset  
*      1: Boot to UART DFU mode  
*      2: Boot to OTA DFU mode
*
* Events generated
*
* gecko_evt_system_boot - Sent after the device has booted in normal mode.
* gecko_evt_dfu_boot - Sent after the device has booted in UART DFU mode.
*
**/

static inline void* gecko_cmd_system_reset(uint8 dfu)
{
    
    gecko_cmd_msg->data.cmd_system_reset.dfu=dfu;
    gecko_cmd_msg->header=(gecko_cmd_system_reset_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command_noresponse(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    return 0;
}

/** 
*
* gecko_cmd_system_get_bt_address
*
* Read the Bluetooth public address used by the device. 
*
*
**/

static inline struct gecko_msg_system_get_bt_address_rsp_t* gecko_cmd_system_get_bt_address()
{
    
    gecko_cmd_msg->header=(gecko_cmd_system_get_bt_address_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_get_bt_address;
}

/** 
*
* gecko_cmd_system_set_bt_address
*
* Deprecated and replaced by system_set_identity_address command.
* 
* Set the Bluetooth public address used by the device. A valid address set with
* this command overrides the default Bluetooth public address programmed at
* production and is effective in the next system reboot. The stack treats
* 00:00:00:00:00:00 and ff:ff:ff:ff:ff:ff as invalid addresses. As a result,
* passing one of them into this command will cause the stack to use the default
* address in the next system reboot. 
*
* @param address   Bluetooth public address in little endian format
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_system_set_bt_address_rsp_t* gecko_cmd_system_set_bt_address(bd_addr address)
{
    
    memcpy(&gecko_cmd_msg->data.cmd_system_set_bt_address.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->header=(gecko_cmd_system_set_bt_address_id+(((6)&0xff)<<8)+(((6)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_set_bt_address;
}

/** 
*
* gecko_cmd_system_set_tx_power
*
* Set the global maximum TX power for Bluetooth. The returned value is the
* selected maximum output power level after applying the RF path compensation.
* If the GATT server contains a TX power service, the TX Power Level attribute
* will be updated accordingly.
* 
* The selected power level may be less than the specified value if the device
* does not meet the power requirements. For Bluetooth connections, the maximum
* TX power is limited to 10 dBm if Adaptive Frequency Hopping (AFH) is not
* enabled.
* 
* By default, the global maximum TX power value is 8 dBm.
* 
* NOTE: Do not use this command while advertising, scanning, or during
* connection. 
*
* @param power   TX power in 0.1 dBm steps. For example, the value of 10 is 1 dBm and 55 is 5.5
*  dBm.
*
**/

static inline struct gecko_msg_system_set_tx_power_rsp_t* gecko_cmd_system_set_tx_power(int16 power)
{
    
    gecko_cmd_msg->data.cmd_system_set_tx_power.power=power;
    gecko_cmd_msg->header=(gecko_cmd_system_set_tx_power_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_set_tx_power;
}

/** 
*
* gecko_cmd_system_get_random_data
*
* Get random data up to 16 bytes. 
*
* @param length   Length of random data. Maximum length is 16 bytes.
*
**/

static inline struct gecko_msg_system_get_random_data_rsp_t* gecko_cmd_system_get_random_data(uint8 length)
{
    
    gecko_cmd_msg->data.cmd_system_get_random_data.length=length;
    gecko_cmd_msg->header=(gecko_cmd_system_get_random_data_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_get_random_data;
}

/** 
*
* gecko_cmd_system_halt
*
* Force radio to idle state and allow device to sleep. Advertising, scanning,
* connections, and software timers are halted by this command. Halted operations
* resume after calling this command with parameter 0. Connections stay alive if
* system is resumed before connection supervision timeout.
* 
* Use this command only for a short time period (a few seconds at maximum).
* Although it halts Bluetooth activity, all tasks and operations still exist
* inside the stack with their own concepts of time. Halting the system for a
* long time period may have negative consequences on stack's internal states.
* 
* NOTE: The software timer is also halted. Hardware interrupts are the only way
* to wake up from energy mode 2 when the system is halted. 
*
* @param halt   Values:
*  
*      1: halt  
*      0: resume
*
**/

static inline struct gecko_msg_system_halt_rsp_t* gecko_cmd_system_halt(uint8 halt)
{
    
    gecko_cmd_msg->data.cmd_system_halt.halt=halt;
    gecko_cmd_msg->header=(gecko_cmd_system_halt_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_halt;
}

/** 
*
* gecko_cmd_system_set_device_name
*
* Set the device name which will be used during the OTA update. The name will be
* stored in the persistent store. If the OTA device name is also set in the
* stack configuration, the name stored in the persistent store is overwritten by
* the name in the stack configuration during the device boot. 
*
* @param type   Device name to set. Values:
*  
*       0: OTA device name
* @param name_len   Array length
* @param name_data   Device name
*
**/

static inline struct gecko_msg_system_set_device_name_rsp_t* gecko_cmd_system_set_device_name(uint8 type,uint8 name_len, const uint8* name_data)
{
    if ((uint16_t)name_len > BGLIB_MSG_MAX_PAYLOAD - 2)
    {
        gecko_rsp_msg->data.rsp_system_set_device_name.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_system_set_device_name;
    }

    
    gecko_cmd_msg->data.cmd_system_set_device_name.type=type;
    gecko_cmd_msg->data.cmd_system_set_device_name.name.len=name_len;
    memcpy(gecko_cmd_msg->data.cmd_system_set_device_name.name.data,name_data,name_len);
    gecko_cmd_msg->header=(gecko_cmd_system_set_device_name_id+(((2+name_len)&0xff)<<8)+(((2+name_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_set_device_name;
}

/** 
*
* gecko_cmd_system_linklayer_configure
*
* Send configuration data to the link layer. This command fine tunes low-level
* Bluetooth operations. 
*
* @param key   Key to configure
* @param data_len   Array length
* @param data_data   Configuration data. Length and contents of the data field depend on the key
*  value used.
*
**/

static inline struct gecko_msg_system_linklayer_configure_rsp_t* gecko_cmd_system_linklayer_configure(uint8 key,uint8 data_len, const uint8* data_data)
{
    if ((uint16_t)data_len > BGLIB_MSG_MAX_PAYLOAD - 2)
    {
        gecko_rsp_msg->data.rsp_system_linklayer_configure.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_system_linklayer_configure;
    }

    
    gecko_cmd_msg->data.cmd_system_linklayer_configure.key=key;
    gecko_cmd_msg->data.cmd_system_linklayer_configure.data.len=data_len;
    memcpy(gecko_cmd_msg->data.cmd_system_linklayer_configure.data.data,data_data,data_len);
    gecko_cmd_msg->header=(gecko_cmd_system_linklayer_configure_id+(((2+data_len)&0xff)<<8)+(((2+data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_linklayer_configure;
}

/** 
*
* gecko_cmd_system_get_counters
*
* Get packet and error counters. Passing a non-zero value also resets counters. 
*
* @param reset   Reset counters if the parameter value is not zero.
*
**/

static inline struct gecko_msg_system_get_counters_rsp_t* gecko_cmd_system_get_counters(uint8 reset)
{
    
    gecko_cmd_msg->data.cmd_system_get_counters.reset=reset;
    gecko_cmd_msg->header=(gecko_cmd_system_get_counters_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_get_counters;
}

/** 
*
* gecko_cmd_system_data_buffer_write
*
* Write data into the system data buffer. Data will be appended to the end of
* existing data. 
*
* @param data_len   Array length
* @param data_data   Data to write
*
**/

static inline struct gecko_msg_system_data_buffer_write_rsp_t* gecko_cmd_system_data_buffer_write(uint8 data_len, const uint8* data_data)
{
    if ((uint16_t)data_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_system_data_buffer_write.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_system_data_buffer_write;
    }

    
    gecko_cmd_msg->data.cmd_system_data_buffer_write.data.len=data_len;
    memcpy(gecko_cmd_msg->data.cmd_system_data_buffer_write.data.data,data_data,data_len);
    gecko_cmd_msg->header=(gecko_cmd_system_data_buffer_write_id+(((1+data_len)&0xff)<<8)+(((1+data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_data_buffer_write;
}

/** 
*
* gecko_cmd_system_set_identity_address
*
* Set the device's Bluetooth identity address. The address can be a public
* device address or a static device address. A valid address set with this
* command will be written into persistent storage using PS keys. The stack
* returns an error if the static device address does not conform to the
* Bluetooth specification.
* 
* The new address will be effective in the next system reboot. The stack will
* use the address in the PS keys when present. Otherwise, it uses the default
* Bluetooth public device address which is programmed at production.
* 
* The stack treats 00:00:00:00:00:00 and ff:ff:ff:ff:ff:ff as invalid addresses.
* Therefore, passing one of them into this command will cause the stack to
* delete the PS keys and use the default address in the next system reboot.
* 
* Note: Because the PS keys are located in flash and flash wearing can occur,
* avoid calling this command regularly. 
*
* @param address   Bluetooth identity address in little endian format
* @param type   Address type
*  
*      0: Public device address  
*      1: Static device address
*
**/

static inline struct gecko_msg_system_set_identity_address_rsp_t* gecko_cmd_system_set_identity_address(bd_addr address,uint8 type)
{
    
    memcpy(&gecko_cmd_msg->data.cmd_system_set_identity_address.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->data.cmd_system_set_identity_address.type=type;
    gecko_cmd_msg->header=(gecko_cmd_system_set_identity_address_id+(((7)&0xff)<<8)+(((7)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_set_identity_address;
}

/** 
*
* gecko_cmd_system_data_buffer_clear
*
* Remove all data from the system data buffer. 
*
*
**/

static inline struct gecko_msg_system_data_buffer_clear_rsp_t* gecko_cmd_system_data_buffer_clear()
{
    
    gecko_cmd_msg->header=(gecko_cmd_system_data_buffer_clear_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_system_data_buffer_clear;
}

/** 
*
* gecko_cmd_le_gap_open
*
* Deprecated and replaced by le_gap_connect command, which allows opening a
* connection with a specified PHY.
* 
* Connect to an advertising device where 1M PHY is the initiating PHY. 
*
* @param address   An address of the device to connect to
* @param address_type   An address type of the device to connect to
*
* Events generated
*
* gecko_evt_le_connection_opened - Triggered after the connection is opened and indicates whether the devices are
*  already bonded and whether the role of the Bluetooth device is Slave or
*  Master.
* gecko_evt_le_connection_parameters - Indicates the connection parameters and security mode of the connection.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_open_rsp_t* gecko_cmd_le_gap_open(bd_addr address,uint8 address_type)
{
    
    memcpy(&gecko_cmd_msg->data.cmd_le_gap_open.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->data.cmd_le_gap_open.address_type=address_type;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_open_id+(((7)&0xff)<<8)+(((7)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_open;
}

/** 
*
* gecko_cmd_le_gap_set_mode
*
* Deprecated. Use le_gap_start_advertising command to enable advertising and
* le_gap_stop_advertising command to disable advertising.
* 
* This command is only effective on the first advertising set (handle value 0).
* Other advertising sets are not affected. 
*
* @param discover   Discoverable mode
* @param connect   Connectable mode
*
* Events generated
*
* gecko_evt_le_gap_adv_timeout - Triggered when the number of advertising events is done and advertising has
*  stopped.
* gecko_evt_le_connection_opened - Triggered when a remote device opens a connection to this advertising device.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_set_mode_rsp_t* gecko_cmd_le_gap_set_mode(uint8 discover,uint8 connect)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_mode.discover=discover;
    gecko_cmd_msg->data.cmd_le_gap_set_mode.connect=connect;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_mode_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_mode;
}

/** 
*
* gecko_cmd_le_gap_discover
*
* Deprecated and replaced by le_gap_start_discovery command. To preserve the
* same functionality when migrating to the new command, use 1M PHY in
* scanning_phy parameter.
* 
* This command can be used to start the GAP discovery procedure to scan for
* advertising devices on 1M PHY. To cancel an ongoing discovery process, use the
* le_gap_end_procedure command. 
*
* @param mode   Bluetooth discovery Mode. For values see link.
*
* Events generated
*
* gecko_evt_le_gap_scan_response - Each time an advertising packet is received, this event is triggered. The
*  packets are not filtered in any way, so multiple events will be received for
*  every advertising device in range.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_discover_rsp_t* gecko_cmd_le_gap_discover(uint8 mode)
{
    
    gecko_cmd_msg->data.cmd_le_gap_discover.mode=mode;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_discover_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_discover;
}

/** 
*
* gecko_cmd_le_gap_end_procedure
*
* End the current GAP discovery procedure (i.e., scanning for advertising
* devices). 
*
*
**/

static inline struct gecko_msg_le_gap_end_procedure_rsp_t* gecko_cmd_le_gap_end_procedure()
{
    
    gecko_cmd_msg->header=(gecko_cmd_le_gap_end_procedure_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_end_procedure;
}

/** 
*
* gecko_cmd_le_gap_set_adv_parameters
*
* Deprecated and replaced by le_gap_set_advertise_timing command to set the
* advertising intervals and le_gap_set_advertise_channel_map command to set the
* channel map.
* 
* This command is only effective on the first advertising set (handle value 0).
* Other advertising sets are not affected. 
*
* @param interval_min   Minimum advertising interval. Value in units of 0.625 ms
*  
*      Range: 0x20 to 0xFFFF  
*      Time range: 20 ms to 40.96 s  
*  
*  Default value: 100 ms
* @param interval_max   Maximum advertising interval. Value in units of 0.625 ms
*  
*      Range: 0x20 to 0xFFFF  
*      Time range: 20 ms to 40.96 s  
*      Note: interval_max should be bigger than interval_min  
*  
*  Default value: 200 ms
* @param channel_map   Advertising channel map, which determines which of the three channels will be
*  used for advertising. This value is given as a bitmask. Values:
*  
*       1: Advertise on CH37  
*       2: Advertise on CH38  
*       3: Advertise on CH37 and CH38  
*       4: Advertise on CH39  
*       5: Advertise on CH37 and CH39  
*       6: Advertise on CH38 and CH39  
*       7: Advertise on all channels  
*  
*  Recommended value: 7
*  
*  Default value: 7
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_set_adv_parameters_rsp_t* gecko_cmd_le_gap_set_adv_parameters(uint16 interval_min,uint16 interval_max,uint8 channel_map)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_adv_parameters.interval_min=interval_min;
    gecko_cmd_msg->data.cmd_le_gap_set_adv_parameters.interval_max=interval_max;
    gecko_cmd_msg->data.cmd_le_gap_set_adv_parameters.channel_map=channel_map;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_adv_parameters_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_adv_parameters;
}

/** 
*
* gecko_cmd_le_gap_set_conn_parameters
*
* Deprecated and replaced by le_gap_set_conn_timing_parameters command for
* setting timing parameters.
* 
* Set the default Bluetooth connection parameters. The configured values are
* valid for all subsequent connections that will be established. To change the
* parameters of an already established connection, use the command
* le_connection_set_parameters. 
*
* @param min_interval   Minimum value for the connection event interval. This must be set less than or
*  equal to the max_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s  
*  
*  Default value: 20 ms
* @param max_interval   Maximum value for the connection event interval. This must be set greater than
*  or equal to the min_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s  
*  
*  Default value: 50 ms
* @param latency   Slave latency, which defines how many connection intervals the slave can skip
*  if it has no data to send
*  
*      Range: 0x0000 to 0x01f4  
*  
*  Default value: 0
* @param timeout   Supervision timeout, which defines the time that the connection is maintained
*  although the devices can't communicate at the currently configured connection
*  intervals.
*  
*      Range: 0x000a to 0x0c80  
*      Time = Value x 10 ms  
*      Time Range: 100 ms to 32 s  
*      The value in milliseconds must be larger than (1 + latency) * max_interval
*      * 2, where max_interval is given in milliseconds  
*  
*  Set the supervision timeout at a value which allows communication attempts
*  over at least a few connection intervals.
*  
*  Default value: 1000 ms
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_set_conn_parameters_rsp_t* gecko_cmd_le_gap_set_conn_parameters(uint16 min_interval,uint16 max_interval,uint16 latency,uint16 timeout)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_conn_parameters.min_interval=min_interval;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_parameters.max_interval=max_interval;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_parameters.latency=latency;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_parameters.timeout=timeout;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_conn_parameters_id+(((8)&0xff)<<8)+(((8)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_conn_parameters;
}

/** 
*
* gecko_cmd_le_gap_set_scan_parameters
*
* Deprecated and replaced by le_gap_set_discovery_timing command to set timing
* parameters, and le_gap_set_discovery_type command for the scan type.
* 
* The parameters set by this command are only effective on the 1M PHY. For Coded
* PHY, use the above replacement command. 
*
* @param scan_interval   Scanner interval is defined as the time interval when the device starts its
*  last scan until it begins the subsequent scan. In other words, it indicates
*  how often to scan
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0004 to 0x4000  
*      Time Range: 2.5 ms to 10.24 s  
*  
*  Default value: 10 ms
*  
*  A variable delay occurs when switching channels at the end of each scanning
*  interval, which is included in the scanning interval time. During the switch
*  time no advertising packets are received by the device. The switch time
*  variation is use case-dependent. For example, if scanning while keeping active
*  connections, the channel switch time might be longer than scanning without any
*  active connections. Increasing the scanning interval reduces the amount of
*  time in which the device can't receive advertising packets because it will
*  switch channels less often.
*  
*  After every scan interval, the scanner changes the frequency at which it
*  operates. It cycles through all three advertising channels in a round robin
*  fashion. According to the specification, all three channels must be used by
*  the scanner.
* @param scan_window   Scan window defines the duration of the scan which must be less than or equal
*  to scan_interval
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0004 to 0x4000  
*      Time Range: 2.5 ms to 10.24 s  
*  
*  Default value: 10 ms Note that packet reception is aborted if it was started
*  before the scan window ends.
* @param active   The scan type. Values:
*  
*       0: Passive scanning  
*       1: Active scanning  
*      In passive scanning mode, the device only listens to advertising packets
*      and does not transmit any packets.  
*      In active scanning mode, the device will send out a scan request packet
*      upon receiving an advertising packet from a remote device and then it will
*      listen to the scan response packet from the device.  
*  
*  Default value: 0
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_set_scan_parameters_rsp_t* gecko_cmd_le_gap_set_scan_parameters(uint16 scan_interval,uint16 scan_window,uint8 active)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_scan_parameters.scan_interval=scan_interval;
    gecko_cmd_msg->data.cmd_le_gap_set_scan_parameters.scan_window=scan_window;
    gecko_cmd_msg->data.cmd_le_gap_set_scan_parameters.active=active;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_scan_parameters_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_scan_parameters;
}

/** 
*
* gecko_cmd_le_gap_set_adv_data
*
* Deprecated. Use le_gap_bt5_set_adv_data command to set advertising data and
* scan response data.
* 
* This command is only effective on the first advertising set (handle value 0).
* Other advertising sets are not affected. 
*
* @param scan_rsp   This value selects if data is intended for advertising packets, scan response
*  packets, or advertising packet in OTA. Values:
*  
*      0: Advertising packets  
*       1: Scan response packets  
*       2: OTA advertising packets  
*       4: OTA scan response packets
* @param adv_data_len   Array length
* @param adv_data_data   Data to be set
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_set_adv_data_rsp_t* gecko_cmd_le_gap_set_adv_data(uint8 scan_rsp,uint8 adv_data_len, const uint8* adv_data_data)
{
    if ((uint16_t)adv_data_len > BGLIB_MSG_MAX_PAYLOAD - 2)
    {
        gecko_rsp_msg->data.rsp_le_gap_set_adv_data.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_le_gap_set_adv_data;
    }

    
    gecko_cmd_msg->data.cmd_le_gap_set_adv_data.scan_rsp=scan_rsp;
    gecko_cmd_msg->data.cmd_le_gap_set_adv_data.adv_data.len=adv_data_len;
    memcpy(gecko_cmd_msg->data.cmd_le_gap_set_adv_data.adv_data.data,adv_data_data,adv_data_len);
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_adv_data_id+(((2+adv_data_len)&0xff)<<8)+(((2+adv_data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_adv_data;
}

/** 
*
* gecko_cmd_le_gap_set_adv_timeout
*
* Deprecated. Use the new command le_gap_set_advertise_timing.
* 
* This command is only effective on the first advertising set (handle value 0).
* Other advertising sets are not affected. 
*
* @param maxevents   If non-zero, indicates the maximum number of advertising events to send before
*  stopping advertiser. Value 0 indicates no maximum number limit.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_set_adv_timeout_rsp_t* gecko_cmd_le_gap_set_adv_timeout(uint8 maxevents)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_adv_timeout.maxevents=maxevents;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_adv_timeout_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_adv_timeout;
}

/** 
*
* gecko_cmd_le_gap_set_conn_phy
*
* Set default preferred and accepted PHYs. PHY settings will be used for all
* subsequent connections. Non-preferred PHY can also be set if the remote device
* does not accept any of the preferred PHYs.
* 
* The parameter accepted_phy is used to specify PHYs that the stack can accept
* in a remotely-initiated PHY update request. A PHY update will not happen if
* none of the accepted PHYs are present in the request.
* 
* NOTE: 2M and Coded PHYs are not supported by all devices. 
*
* @param preferred_phy   Preferred PHYs. This parameter is a bitfield and multiple PHYs can be set.
*  
*      0x01: 1M PHY  
*      0x02: 2M PHY  
*      0x04: Coded PHY  
*      0xff: Any PHYs  
*  
*  Default: 0xff (no preference)
* @param accepted_phy   Accepted PHYs in remotely-initiated PHY update request. This parameter is a
*  bitfield and multiple PHYs can be set.
*  
*      0x01: 1M PHY  
*      0x02: 2M PHY  
*      0x04: Coded PHY  
*      0xff: Any PHYs  
*  
*  Default: 0xff (all PHYs accepted)
*
**/

static inline struct gecko_msg_le_gap_set_conn_phy_rsp_t* gecko_cmd_le_gap_set_conn_phy(uint8 preferred_phy,uint8 accepted_phy)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_conn_phy.preferred_phy=preferred_phy;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_phy.accepted_phy=accepted_phy;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_conn_phy_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_conn_phy;
}

/** 
*
* gecko_cmd_le_gap_bt5_set_mode
*
* Deprecated and replaced by le_gap_start_advertising command to start
* advertising, and le_gap_stop_advertising command to stop advertising.
* le_gap_set_advertise_timing command can be used for setting the maxevents and
* command le_gap_set_advertise_configuration can be used for setting address
* types. 
*
* @param handle   Advertising set handle
* @param discover   Discoverable mode
* @param connect   Connectable mode
* @param maxevents   If non-zero, indicates the maximum number of advertising events to send before
*  stopping the advertiser. Value 0 indicates no maximum number limit.
* @param address_type   Address type to use for packets
*
* Events generated
*
* gecko_evt_le_gap_adv_timeout - Triggered when the advertising events set by this command are complete and
*  advertising is stopped on the given advertising set.
* gecko_evt_le_connection_opened - Triggered when a remote device opens a connection to the advertiser on the
*  specified advertising set.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_bt5_set_mode_rsp_t* gecko_cmd_le_gap_bt5_set_mode(uint8 handle,uint8 discover,uint8 connect,uint16 maxevents,uint8 address_type)
{
    
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_mode.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_mode.discover=discover;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_mode.connect=connect;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_mode.maxevents=maxevents;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_mode.address_type=address_type;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_bt5_set_mode_id+(((6)&0xff)<<8)+(((6)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_bt5_set_mode;
}

/** 
*
* gecko_cmd_le_gap_bt5_set_adv_parameters
*
* Deprecated and replaced by le_gap_set_advertise_timing command to set the
* advertising intervals, le_gap_set_advertise_channel_map command to set the
* channel map, and le_gap_set_advertise_report_scan_request command to enable
* and disable scan request notifications. 
*
* @param handle   Advertising set handle
* @param interval_min   Minimum advertising interval. Value in units of 0.625 ms
*  
*      Range: 0x20 to 0xFFFF  
*      Time range: 20 ms to 40.96 s  
*  
*  Default value: 100 ms
* @param interval_max   Maximum advertising interval. Value in units of 0.625 ms
*  
*      Range: 0x20 to 0xFFFF  
*      Time range: 20 ms to 40.96 s  
*      Note: interval_max should be bigger than interval_min  
*  
*  Default value: 200 ms
* @param channel_map   Advertising channel map, which determines which of the three channels will be
*  used for advertising. This value is given as a bitmask. Values:
*  
*       1: Advertise on CH37  
*       2: Advertise on CH38  
*       3: Advertise on CH37 and CH38  
*       4: Advertise on CH39  
*       5: Advertise on CH37 and CH39  
*       6: Advertise on CH38 and CH39  
*       7: Advertise on all channels  
*  
*  Recommended value: 7
*  
*  Default value: 7
* @param report_scan   If non-zero, enables scan request notification, and scan requests will be
*  reported as events.
*  
*  Default value: 0
*
* Events generated
*
* gecko_evt_le_gap_scan_request - Triggered when a scan request is received during advertising if the scan
*  request notification is enabled by this command.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_gap_bt5_set_adv_parameters_rsp_t* gecko_cmd_le_gap_bt5_set_adv_parameters(uint8 handle,uint16 interval_min,uint16 interval_max,uint8 channel_map,uint8 report_scan)
{
    
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_parameters.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_parameters.interval_min=interval_min;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_parameters.interval_max=interval_max;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_parameters.channel_map=channel_map;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_parameters.report_scan=report_scan;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_bt5_set_adv_parameters_id+(((7)&0xff)<<8)+(((7)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_bt5_set_adv_parameters;
}

/** 
*
* gecko_cmd_le_gap_bt5_set_adv_data
*
* Set user-defined data in advertising packets, scan response packets, or
* periodic advertising packets. Maximum 31 bytes of data can be set for legacy
* advertising. Maximum 191 bytes of data can be set for connectable extended
* advertising. Maximum 253 bytes of data can be set for periodic and non-
* connectable extended advertising. For setting longer advertising data, use
* command le_gap_set_long_advertising_data.
* 
* If advertising mode is currently enabled, the new advertising data will be
* used immediately. Advertising mode can be enabled using command
* le_gap_start_advertising. Periodic advertising mode can be enabled using
* command le_gap_start_periodic_advertising.
* 
* The invalid parameter error will be returned in the following situations:
* 
*      Data length is more than 31 bytes but the advertiser can only advertise
*     using legacy advertising PDUs.  
*      Data is too long to fit into a single advertisement.  
*      Set data of the advertising data packet when the advertiser is advertising
*     in scannable mode using extended advertising PDUs.  
*      Set data of the scan response data packet when the advertiser is
*     advertising in connectable mode using extended advertising PDUs.  
* 
* Note that the user-defined data may be overwritten by the system when the
* advertising is later enabled in a discoverable mode other than user_data. 
*
* @param handle   Advertising set handle
* @param scan_rsp   This value selects whether data is intended for advertising packets, scan
*  response packets, periodic advertising packets, or advertising packets in OTA.
*  Values are as follows:
*  
*      0: Advertising packets  
*      1: Scan response packets  
*      2: OTA advertising packets  
*      4: OTA scan response packets  
*      8: Periodic advertising packets
* @param adv_data_len   Array length
* @param adv_data_data   Data to be set
*
**/

static inline struct gecko_msg_le_gap_bt5_set_adv_data_rsp_t* gecko_cmd_le_gap_bt5_set_adv_data(uint8 handle,uint8 scan_rsp,uint8 adv_data_len, const uint8* adv_data_data)
{
    if ((uint16_t)adv_data_len > BGLIB_MSG_MAX_PAYLOAD - 3)
    {
        gecko_rsp_msg->data.rsp_le_gap_bt5_set_adv_data.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_le_gap_bt5_set_adv_data;
    }

    
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_data.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_data.scan_rsp=scan_rsp;
    gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_data.adv_data.len=adv_data_len;
    memcpy(gecko_cmd_msg->data.cmd_le_gap_bt5_set_adv_data.adv_data.data,adv_data_data,adv_data_len);
    gecko_cmd_msg->header=(gecko_cmd_le_gap_bt5_set_adv_data_id+(((3+adv_data_len)&0xff)<<8)+(((3+adv_data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_bt5_set_adv_data;
}

/** 
*
* gecko_cmd_le_gap_set_privacy_mode
*
* Enable or disable the privacy feature on all GAP roles. New privacy mode will
* take effect for advertising next time advertising is enabled, for scanning
* next time scanning is enabled, and for initiating on the next open connection
* command. When privacy is enabled and the device is advertising or scanning,
* the stack will maintain a periodic timer with the specified time interval as a
* timeout value. At each timeout, the stack will generate a new private
* resolvable address and use it in advertising data packets and scanning
* requests.
* 
* By default, privacy feature is disabled. 
*
* @param privacy   Values:
*  
*      0: Disable privacy  
*       1: Enable privacy
* @param interval   The minimum time interval between a private address change. This parameter is
*  ignored if this command is issued to disable privacy mode. Values:
*  
*      0: Use default interval, 15 minutes  
*       others: The time interval in minutes
*
**/

static inline struct gecko_msg_le_gap_set_privacy_mode_rsp_t* gecko_cmd_le_gap_set_privacy_mode(uint8 privacy,uint8 interval)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_privacy_mode.privacy=privacy;
    gecko_cmd_msg->data.cmd_le_gap_set_privacy_mode.interval=interval;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_privacy_mode_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_privacy_mode;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_timing
*
* Set the advertising timing parameters of the given advertising set. This
* setting will take effect next time that advertising is enabled. 
*
* @param handle   Advertising set handle
* @param interval_min   Minimum advertising interval. Value in units of 0.625 ms
*  
*      Range: 0x20 to 0xFFFF  
*      Time range: 20 ms to 40.96 s  
*  
*  Default value: 100 ms
* @param interval_max   Maximum advertising interval. Value in units of 0.625 ms
*  
*      Range: 0x20 to 0xFFFF  
*      Time range: 20 ms to 40.96 s  
*      Note: interval_max should be bigger than interval_min  
*  
*  Default value: 200 ms
* @param duration   Advertising duration for this advertising set. Value 0 indicates no
*  advertising duration limit and advertising continues until it is disabled. A
*  non-zero value sets the duration in units of 10 ms. The duration begins at the
*  start of the first advertising event of this advertising set.
*  
*      Range: 0x0001 to 0xFFFF  
*      Time range: 10 ms to 655.35 s  
*  
*  Default value: 0
* @param maxevents   If non-zero, indicates the maximum number of advertising events to send before
*  the advertiser is stopped. Value 0 indicates no maximum number limit.
*  
*  Default value: 0
*
**/

static inline struct gecko_msg_le_gap_set_advertise_timing_rsp_t* gecko_cmd_le_gap_set_advertise_timing(uint8 handle,uint32 interval_min,uint32 interval_max,uint16 duration,uint8 maxevents)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_timing.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_timing.interval_min=interval_min;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_timing.interval_max=interval_max;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_timing.duration=duration;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_timing.maxevents=maxevents;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_timing_id+(((12)&0xff)<<8)+(((12)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_timing;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_channel_map
*
* Set the primary advertising channel map of the given advertising set. This
* setting will take effect next time that advertising is enabled. 
*
* @param handle   Advertising set handle
* @param channel_map   Advertising channel map which determines which of the three channels will be
*  used for advertising. This value is given as a bitmask. Values:
*  
*       1: Advertise on CH37  
*       2: Advertise on CH38  
*       3: Advertise on CH37 and CH38  
*       4: Advertise on CH39  
*       5: Advertise on CH37 and CH39  
*       6: Advertise on CH38 and CH39  
*       7: Advertise on all channels  
*  
*  Recommended value: 7
*  
*  Default value: 7
*
**/

static inline struct gecko_msg_le_gap_set_advertise_channel_map_rsp_t* gecko_cmd_le_gap_set_advertise_channel_map(uint8 handle,uint8 channel_map)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_channel_map.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_channel_map.channel_map=channel_map;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_channel_map_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_channel_map;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_report_scan_request
*
* Enable or disable the scan request notification of a given advertising set.
* This setting will take effect next time that advertising is enabled. 
*
* @param handle   Advertising set handle
* @param report_scan_req   If non-zero, enables scan request notification and scan requests will be
*  reported as events.
*  
*  Default value: 0
*
* Events generated
*
* gecko_evt_le_gap_scan_request - Triggered when a scan request is received during advertising if the scan
*  request notification is enabled by this command.
*
**/

static inline struct gecko_msg_le_gap_set_advertise_report_scan_request_rsp_t* gecko_cmd_le_gap_set_advertise_report_scan_request(uint8 handle,uint8 report_scan_req)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_report_scan_request.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_report_scan_request.report_scan_req=report_scan_req;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_report_scan_request_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_report_scan_request;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_phy
*
* Set advertising PHYs of the given advertising set. This setting will take
* effect next time that advertising is enabled. The invalid parameter error is
* returned if a PHY value is invalid or the device does not support a given PHY. 
*
* @param handle   Advertising set handle
* @param primary_phy   The PHY on which the advertising packets are transmitted on the primary
*  advertising channel. If legacy advertising PDUs are used, 1M PHY must be used.
*  
*  Values:
*  
*       1: Advertising PHY is 1M PHY  
*       4: Advertising PHY is Coded PHY  
*  
*  Default: 1
* @param secondary_phy   The PHY on which the advertising packets are transmitted on the secondary
*  advertising channel.
*  
*  Values:
*  
*       1: Advertising PHY is 1M PHY  
*       2: Advertising PHY is 2M PHY  
*       4: Advertising PHY is Coded PHY  
*  
*  Default: 1
*
**/

static inline struct gecko_msg_le_gap_set_advertise_phy_rsp_t* gecko_cmd_le_gap_set_advertise_phy(uint8 handle,uint8 primary_phy,uint8 secondary_phy)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_phy.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_phy.primary_phy=primary_phy;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_phy.secondary_phy=secondary_phy;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_phy_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_phy;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_configuration
*
* Enable advertising configuration flags on the given advertising set. The
* configuration change will take effect next time that advertising is enabled.
* 
* These configuration flags can be disabled using
* le_gap_clear_advertise_configuration command. 
*
* @param handle   Advertising set handle
* @param configurations   Advertising configuration flags to enable. This value can be a bitmask of
*  multiple flags. Flags:
*  
*      1 (Bit 0): Use legacy advertising PDUs.  
*      2 (Bit 1): Omit advertiser's address from all PDUs (anonymous
*      advertising). This flag is effective only in extended advertising.  
*       4 (Bit 2): Use le_gap_non_resolvable address type. Advertising must be in
*      non-connectable mode if this configuration is enabled.  
*       8 (Bit 3): Include TX power in advertising packets. This flag is
*      effective only in extended advertising.  
*  
*  Default value: 1
*
**/

static inline struct gecko_msg_le_gap_set_advertise_configuration_rsp_t* gecko_cmd_le_gap_set_advertise_configuration(uint8 handle,uint32 configurations)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_configuration.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_configuration.configurations=configurations;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_configuration_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_configuration;
}

/** 
*
* gecko_cmd_le_gap_clear_advertise_configuration
*
* Disable advertising configuration flags on the given advertising set. The
* configuration change will take effect next time that advertising is enabled.
* 
* These configuration flags can be enabled using
* le_gap_set_advertise_configuration command. 
*
* @param handle   Advertising set handle
* @param configurations   Advertising configuration flags to disable. This value can be a bitmask of
*  multiple flags. See le_gap_set_advertise_configuration for possible flags.
*
**/

static inline struct gecko_msg_le_gap_clear_advertise_configuration_rsp_t* gecko_cmd_le_gap_clear_advertise_configuration(uint8 handle,uint32 configurations)
{
    
    gecko_cmd_msg->data.cmd_le_gap_clear_advertise_configuration.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_clear_advertise_configuration.configurations=configurations;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_clear_advertise_configuration_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_clear_advertise_configuration;
}

/** 
*
* gecko_cmd_le_gap_start_advertising
*
* Start advertising of a given advertising set with specified discoverable and
* connectable modes.
* 
* The number of concurrent advertising is limited by MAX_ADVERTISERS
* configuration.
* 
* The number of concurrent connectable advertising is also limited by
* MAX_CONNECTIONS configuration. For example, only one connectable advertising
* can be enabled if the device has (MAX_CONNECTIONS - 1) connections when this
* command is called. The limitation does not apply to non-connectable
* advertising.
* 
* The default advertising configuration in the stack is set to using legacy
* advertising PDUs on 1M PHY. The stack will automatically select extended
* advertising PDUs if either of the following has occurred with the default
* configuration:
* 
*   1. The connectable mode is set to le_gap_connectable_non_scannable.
*   2. The primary advertising PHY is set to Coded PHY by the command le_gap_set_advertise_phy.
*   3. The user advertising data length is more than 31 bytes.
*   4. Periodic advertising is enabled.
* 
* If currently set parameters can't be used, an error is returned. Specifically,
* this command fails with the connection limit exceeded error if it causes the
* number of connections exceeding the configured MAX_CONNECTIONS value. It fails
* with the invalid parameter error if one of the following use cases occurs:
* 
*   1. Non-resolvable random address is used but the connectable mode is le_gap_connectable_scannable or le_gap_connectable_non_scannable.
*   2. le_gap_connectable_non_scannable is the connectable mode but using legacy advertising PDUs has been explicitly enabled with command le_gap_set_advertise_configuration.
*   3. Coded PHY is the primary advertising PHY but using legacy advertising PDUs has been explicitly enabled with command le_gap_set_advertise_configuration.
*   4. le_gap_connectable_scannable is the connectable mode but using extended advertising PDUs has been explicitly enabled or the primary advertising PHY is set to Coded PHY.
* 
* If advertising is enabled in user_data mode, use le_gap_bt5_set_adv_data to
* set advertising and scan response data before issuing this command. When
* advertising is enabled in modes other than user_data, advertising and scan
* response data is generated by the stack using the following procedure:
* 
*   1. Add a flags field to advertising data.
*   2. Add a TX power level field to advertising data if the TX power service exists in the local GATT database.
*   3. Add a slave connection interval range field to advertising data if the GAP peripheral preferred connection parameters characteristic exists in the local GATT database.
*   4. Add a list of 16-bit service UUIDs to advertising data if there are one or more 16-bit service UUIDs to advertise. The list is complete if all advertised 16-bit UUIDs are in advertising data. Otherwise, the list is incomplete.
*   5. Add a list of 128-bit service UUIDs to advertising data if there are one or more 128-bit service UUIDs to advertise and there is still free space for this field. The list is complete if all advertised 128-bit UUIDs are in advertising data. Otherwise, the list is incomplete. Note that an advertising data packet can contain at most one 128-bit service UUID.
*   6. Try to add the full local name to advertising data if the device is not in privacy mode. If the full local name does not fit into the remaining free space, the advertised name is a shortened version by cutting off the end if the free space has at least 6 bytes. Otherwise, the local name is added to scan response data.
* 
* Event le_connection_opened will be received when a remote device opens a
* connection to the advertiser on this advertising set and also advertising on
* the given set stops.
* 
* Event le_gap_adv_timeout will be received when the number of advertising
* events set by le_gap_set_advertise_timing command is done and advertising with
* the current set has stopped. 
*
* @param handle   Advertising set handle
* @param discover   Discoverable mode
* @param connect   Connectable mode
*
* Events generated
*
* gecko_evt_le_gap_adv_timeout - Triggered when the number of advertising events set by
*  le_gap_set_advertise_timing command is done and advertising has stopped on the
*  given advertising set.
* gecko_evt_le_connection_opened - Triggered when a remote device opens a connection to the advertiser on the
*  specified advertising set and also advertising with the current set stops.
*
**/

static inline struct gecko_msg_le_gap_start_advertising_rsp_t* gecko_cmd_le_gap_start_advertising(uint8 handle,uint8 discover,uint8 connect)
{
    
    gecko_cmd_msg->data.cmd_le_gap_start_advertising.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_start_advertising.discover=discover;
    gecko_cmd_msg->data.cmd_le_gap_start_advertising.connect=connect;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_start_advertising_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_start_advertising;
}

/** 
*
* gecko_cmd_le_gap_stop_advertising
*
* Stop the advertising of the given advertising set. 
*
* @param handle   Advertising set handle
*
**/

static inline struct gecko_msg_le_gap_stop_advertising_rsp_t* gecko_cmd_le_gap_stop_advertising(uint8 handle)
{
    
    gecko_cmd_msg->data.cmd_le_gap_stop_advertising.handle=handle;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_stop_advertising_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_stop_advertising;
}

/** 
*
* gecko_cmd_le_gap_set_discovery_timing
*
* Set the timing parameters of scanning on the specified PHYs. If the device is
* currently scanning for advertising devices on PHYs, new parameters will take
* effect when scanning is restarted. 
*
* @param phys   PHYs for which the parameters are set.
*  
*      1: 1M PHY  
*       4: Coded PHY  
*       5: 1M PHY and Coded PHY
* @param scan_interval   Scan interval is defined as the time interval when the device starts its last
*  scan until it begins the subsequent scan. In other words, how often to scan
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0004 to 0xFFFF  
*      Time Range: 2.5 ms to 40.96 s  
*  
*  Default value: 10 ms
*  
*  A variable delay occurs when switching channels at the end of each scanning
*  interval, which is included in the scanning interval time. During the switch
*  time, advertising packets are not received by the device. The switch time
*  variation is use case-dependent. For example, if scanning while keeping active
*  connections, the channel switch time might be longer than when scanning
*  without any active connections. Increasing the scanning interval reduces the
*  amount of time in which the device can't receive advertising packets because
*  it switches channels less often.
*  
*  After every scan interval, the scanner changes the frequency at which it
*  operates. It cycles through all three advertising channels in a round robin
*  fashion. According to the specification, all three channels must be used by a
*  scanner.
* @param scan_window   Scan window defines the duration of the scan which must be less than or equal
*  to the scan_interval
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0004 to 0xFFFF  
*      Time Range: 2.5 ms to 40.96 s  
*  
*  Default value: 10 ms Note that the packet reception is aborted if it's started
*  just before the scan window ends.
*
**/

static inline struct gecko_msg_le_gap_set_discovery_timing_rsp_t* gecko_cmd_le_gap_set_discovery_timing(uint8 phys,uint16 scan_interval,uint16 scan_window)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_discovery_timing.phys=phys;
    gecko_cmd_msg->data.cmd_le_gap_set_discovery_timing.scan_interval=scan_interval;
    gecko_cmd_msg->data.cmd_le_gap_set_discovery_timing.scan_window=scan_window;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_discovery_timing_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_discovery_timing;
}

/** 
*
* gecko_cmd_le_gap_set_discovery_type
*
* Set the scan type on the specified PHYs. If the device is currently scanning
* for advertising devices on PHYs, new parameters will take effect when scanning
* is restarted. 
*
* @param phys   PHYs for which the parameters are set.
*  
*      1: 1M PHY  
*       4: Coded PHY  
*       5: 1M PHY and Coded PHY
* @param scan_type   Scan type. Values:
*  
*       0: Passive scanning  
*       1: Active scanning  
*      In passive scanning mode, the device only listens to advertising packets
*      and does not transmit packets.  
*      In active scanning mode, the device sends out a scan request packet upon
*      receiving an advertising packet from a remote device. Then, it listens to
*      the scan response packet from the remote device.  
*  
*  Default value: 0
*
**/

static inline struct gecko_msg_le_gap_set_discovery_type_rsp_t* gecko_cmd_le_gap_set_discovery_type(uint8 phys,uint8 scan_type)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_discovery_type.phys=phys;
    gecko_cmd_msg->data.cmd_le_gap_set_discovery_type.scan_type=scan_type;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_discovery_type_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_discovery_type;
}

/** 
*
* gecko_cmd_le_gap_start_discovery
*
* Start the GAP discovery procedure to scan for advertising devices on the
* specified scanning PHY or to perform a device discovery. To cancel an ongoing
* discovery process use the le_gap_end_procedure command.
* 
* The invalid parameter error will be returned if the scanning PHY value is
* invalid or the device does not support the PHY. 
*
* @param scanning_phy   The scanning PHY. Value:
*  
*      1: 1M PHY  
*       4: Coded PHY
* @param mode   Bluetooth discovery Mode. For values see link.
*
* Events generated
*
* gecko_evt_le_gap_scan_response - This event is triggered each time an advertising packet is received. Packets
*  are not filtered in any way, so multiple events will be received for every
*  advertising device in range.
* gecko_evt_le_gap_extended_scan_response - Each time an advertising packet is received and extended scan response is
*  enabled (by le_gap_set_discovery_extended_scan_response), this event is
*  triggered. The packets are not filtered in any way. As a result, multiple
*  events will be received for every advertising device in range.
*
**/

static inline struct gecko_msg_le_gap_start_discovery_rsp_t* gecko_cmd_le_gap_start_discovery(uint8 scanning_phy,uint8 mode)
{
    
    gecko_cmd_msg->data.cmd_le_gap_start_discovery.scanning_phy=scanning_phy;
    gecko_cmd_msg->data.cmd_le_gap_start_discovery.mode=mode;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_start_discovery_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_start_discovery;
}

/** 
*
* gecko_cmd_le_gap_set_data_channel_classification
*
* Specify a channel classification for data channels. This classification
* persists until overwritten with a subsequent command or until the system is
* reset. 
*
* @param channel_map_len   Array length
* @param channel_map_data   This parameter is 5 bytes and contains 37 1-bit fields.  
*  The nth such field (in the range 0 to 36) contains the value for the link
*  layer channel index n.  
*  
*      0: Channel n is bad.  
*       1: Channel n is unknown.  
*  
*  The rest of most significant bits are reserved for future use and shall be set
*  to 0.  
*  At least two channels shall be marked as unknown.
*
**/

static inline struct gecko_msg_le_gap_set_data_channel_classification_rsp_t* gecko_cmd_le_gap_set_data_channel_classification(uint8 channel_map_len, const uint8* channel_map_data)
{
    if ((uint16_t)channel_map_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_le_gap_set_data_channel_classification.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_le_gap_set_data_channel_classification;
    }

    
    gecko_cmd_msg->data.cmd_le_gap_set_data_channel_classification.channel_map.len=channel_map_len;
    memcpy(gecko_cmd_msg->data.cmd_le_gap_set_data_channel_classification.channel_map.data,channel_map_data,channel_map_len);
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_data_channel_classification_id+(((1+channel_map_len)&0xff)<<8)+(((1+channel_map_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_data_channel_classification;
}

/** 
*
* gecko_cmd_le_gap_connect
*
* Connect to an advertising device with the specified initiating PHY on which
* connectable advertisements on primary advertising channels are received. The
* Bluetooth stack will enter a state where it continuously scans for the
* connectable advertising packets from the remote device, which matches the
* Bluetooth address given as a parameter. Scan parameters set in
* le_gap_set_discovery_timing are used in this operation. Upon receiving the
* advertising packet, the module will send a connection request packet to the
* target device to initiate a Bluetooth connection. To cancel an ongoing
* connection process, use the le_connection_close command with the handle
* received in response from this command.
* 
* A connection is opened in no-security mode. If the GATT client needs to read
* or write the attributes on GATT server requiring encryption or authentication,
* it must first encrypt the connection using an appropriate authentication
* method.
* 
* If a connection can't be established (for example, the remote device has gone
* out of range, has entered into deep sleep, or is not advertising), the stack
* will try to connect forever. In this case, the application will not get an
* event related to the connection request. To recover from this situation, the
* application can implement a timeout and call le_connection_close to cancel the
* connection request.
* 
* This command fails with the connection limit exceeded error if the number of
* connections attempted exceeds the configured MAX_CONNECTIONS value.
* 
* This command fails with the invalid parameter error if the initiating PHY
* value is invalid or the device does not support PHY.
* 
* Later calls of this command have to wait for the ongoing command to complete.
* A received event le_connection_opened indicates that the connection opened
* successfully and a received event le_connection_closed indicates that
* connection failures have occurred. 
*
* @param address   Address of the device to connect to
* @param address_type   Address type of the device to connect to
* @param initiating_phy   The initiating PHY. Value:
*  
*      1: 1M PHY  
*       4: Coded PHY
*
* Events generated
*
* gecko_evt_le_connection_opened - This event is triggered after the connection is opened and indicates whether
*  the devices are already bonded and whether the role of the Bluetooth device is
*  Slave or Master.
* gecko_evt_le_connection_parameters - This event indicates the connection parameters and security mode of the
*  connection.
*
**/

static inline struct gecko_msg_le_gap_connect_rsp_t* gecko_cmd_le_gap_connect(bd_addr address,uint8 address_type,uint8 initiating_phy)
{
    
    memcpy(&gecko_cmd_msg->data.cmd_le_gap_connect.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->data.cmd_le_gap_connect.address_type=address_type;
    gecko_cmd_msg->data.cmd_le_gap_connect.initiating_phy=initiating_phy;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_connect_id+(((8)&0xff)<<8)+(((8)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_connect;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_tx_power
*
* Limit the maximum advertising TX power on the given advertising set. If the
* value goes over the global value that was set using system_set_tx_power
* command, the global value will be the maximum limit. The maximum TX power of
* legacy advertising is further constrained to be less than +10 dBm. Extended
* advertising TX power can be +10 dBm and over if Adaptive Frequency Hopping is
* enabled.
* 
* This setting will take effect next time advertising is enabled.
* 
* By default, maximum advertising TX power is limited by the global value. 
*
* @param handle   Advertising set handle
* @param power   TX power in 0.1 dBm steps. For example, the value of 10 is 1 dBm and 55 is 5.5
*  dBm.
*
**/

static inline struct gecko_msg_le_gap_set_advertise_tx_power_rsp_t* gecko_cmd_le_gap_set_advertise_tx_power(uint8 handle,int16 power)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_tx_power.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_tx_power.power=power;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_tx_power_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_tx_power;
}

/** 
*
* gecko_cmd_le_gap_set_discovery_extended_scan_response
*
* Enable or disable the extended scan response event. When the extended scan
* response event is enabled, it replaces le_gap_scan_response, that is, the
* stack will generate either le_gap_extended_scan_response or
* le_gap_scan_response, but not both. 
*
* @param enable   Values:
*  
*      0: Disable extended scan response event  
*       1: Enable extended scan response event
*
**/

static inline struct gecko_msg_le_gap_set_discovery_extended_scan_response_rsp_t* gecko_cmd_le_gap_set_discovery_extended_scan_response(uint8 enable)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_discovery_extended_scan_response.enable=enable;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_discovery_extended_scan_response_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_discovery_extended_scan_response;
}

/** 
*
* gecko_cmd_le_gap_start_periodic_advertising
*
* Start periodic advertising on the given advertising set. The stack enables the
* advertising set automatically if the set was not enabled and the set can
* advertise using extended advertising PDUs beside the syncInfo (which is needed
* for the periodic advertising).
* 
* The invalid parameter error is returned if the application has configured
* legacy advertising PDUs or anonymous advertising, or the advertising set is
* enabled using legacy advertising PDUs. 
*
* @param handle   Advertising set handle
* @param interval_min   Minimum periodic advertising interval. Value in units of 1.25 ms
*  
*      Range: 0x06 to 0xFFFF  
*      Time range: 7.5 ms to 81.92 s  
*  
*  Default value: 100 ms
* @param interval_max   Maximum periodic advertising interval. Value in units of 1.25 ms
*  
*      Range: 0x06 to 0xFFFF  
*      Time range: 7.5 ms to 81.92 s  
*      Note: interval_max should be bigger than interval_min  
*  
*  Default value: 200 ms
* @param flags   Periodic advertising configurations. Bitmask of the following:
*  
*      Bit 0: Include TX power in advertising PDU
*
**/

static inline struct gecko_msg_le_gap_start_periodic_advertising_rsp_t* gecko_cmd_le_gap_start_periodic_advertising(uint8 handle,uint16 interval_min,uint16 interval_max,uint32 flags)
{
    
    gecko_cmd_msg->data.cmd_le_gap_start_periodic_advertising.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_start_periodic_advertising.interval_min=interval_min;
    gecko_cmd_msg->data.cmd_le_gap_start_periodic_advertising.interval_max=interval_max;
    gecko_cmd_msg->data.cmd_le_gap_start_periodic_advertising.flags=flags;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_start_periodic_advertising_id+(((9)&0xff)<<8)+(((9)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_start_periodic_advertising;
}

/** 
*
* gecko_cmd_le_gap_stop_periodic_advertising
*
* Stop the periodic advertising on the given advertising set.
* 
* This command does not affect the enable state of the advertising set, i.e.,
* legacy or extended advertising is not stopped. 
*
* @param handle   Advertising set handle
*
**/

static inline struct gecko_msg_le_gap_stop_periodic_advertising_rsp_t* gecko_cmd_le_gap_stop_periodic_advertising(uint8 handle)
{
    
    gecko_cmd_msg->data.cmd_le_gap_stop_periodic_advertising.handle=handle;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_stop_periodic_advertising_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_stop_periodic_advertising;
}

/** 
*
* gecko_cmd_le_gap_set_long_advertising_data
*
* Set advertising data for a specified packet type and advertising set. Data
* currently in the system data buffer will be extracted as the advertising data.
* The buffer will be emptied after this command regardless of the completion
* status.
* 
* Prior to calling this command, add data to the buffer with one or multiple
* calls of system_data_buffer_write.
* 
* Maximum 31 bytes of data can be set for legacy advertising. Maximum 191 bytes
* of data can be set for connectable extended advertising. Maximum 1650 bytes of
* data can be set for periodic and non-connectable extended advertising, but
* advertising parameters may limit the amount of data that can be sent in a
* single advertisement.
* 
* See le_gap_bt5_set_adv_data for more details on advertising data. 
*
* @param handle   Advertising set handle
* @param packet_type   This value selects whether data is intended for advertising packets, scan
*  response packets, or periodic advertising packets. Values:
*  
*      0: Advertising packets  
*      1: Scan response packets  
*      2: OTA advertising packets  
*      4: OTA scan response packets  
*      8: Periodic advertising packets
*
**/

static inline struct gecko_msg_le_gap_set_long_advertising_data_rsp_t* gecko_cmd_le_gap_set_long_advertising_data(uint8 handle,uint8 packet_type)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_long_advertising_data.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_long_advertising_data.packet_type=packet_type;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_long_advertising_data_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_long_advertising_data;
}

/** 
*
* gecko_cmd_le_gap_enable_whitelisting
*
* Enable or disable whitelisting. The setting will be effective the next time
* that scanning is enabled. To add devices to the whitelist, either bond with
* the device or add it manually with sm_add_to_whitelist. 
*
* @param enable   1 enable, 0 disable whitelisting.
*
**/

static inline struct gecko_msg_le_gap_enable_whitelisting_rsp_t* gecko_cmd_le_gap_enable_whitelisting(uint8 enable)
{
    
    gecko_cmd_msg->data.cmd_le_gap_enable_whitelisting.enable=enable;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_enable_whitelisting_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_enable_whitelisting;
}

/** 
*
* gecko_cmd_le_gap_set_conn_timing_parameters
*
* Set the default Bluetooth connection parameters. The configured values are
* valid for all subsequent connections that will be established. To change
* parameters of an already established connection, use the command
* le_connection_set_timing_parameters. 
*
* @param min_interval   Minimum value for the connection event interval. This must be set less than or
*  equal to max_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s  
*  
*  Default value: 20 ms
* @param max_interval   Maximum value for the connection event interval. This must be set greater than
*  or equal to min_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s  
*  
*  Default value: 50 ms
* @param latency   Slave latency, which defines how many connection intervals the slave can skip
*  if it has no data to send
*  
*      Range: 0x0000 to 0x01f4  
*  
*  Default value: 0
* @param timeout   Supervision timeout, which defines the time that the connection is maintained
*  although the devices can't communicate at the currently configured connection
*  intervals.
*  
*      Range: 0x000a to 0x0c80  
*      Time = Value x 10 ms  
*      Time Range: 100 ms to 32 s  
*      The value in milliseconds must be larger than (1 + latency) * max_interval
*      * 2, where max_interval is given in milliseconds  
*  
*  Set the supervision timeout at a value which allows communication attempts
*  over at least a few connection intervals.
*  
*  Default value: 1000 ms
* @param min_ce_length   Minimum value for the connection event length. This must be set be less than
*  or equal to max_ce_length.
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0000 to 0xffff  
*  
*  Default value: 0x0000
*  
*  Value is not currently used and is reserved for future. Set to 0.
* @param max_ce_length   Maximum value for the connection event length. This must be set greater than
*  or equal to min_ce_length.
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0000 to 0xffff  
*  
*  Default value: 0xffff
*
**/

static inline struct gecko_msg_le_gap_set_conn_timing_parameters_rsp_t* gecko_cmd_le_gap_set_conn_timing_parameters(uint16 min_interval,uint16 max_interval,uint16 latency,uint16 timeout,uint16 min_ce_length,uint16 max_ce_length)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_conn_timing_parameters.min_interval=min_interval;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_timing_parameters.max_interval=max_interval;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_timing_parameters.latency=latency;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_timing_parameters.timeout=timeout;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_timing_parameters.min_ce_length=min_ce_length;
    gecko_cmd_msg->data.cmd_le_gap_set_conn_timing_parameters.max_ce_length=max_ce_length;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_conn_timing_parameters_id+(((12)&0xff)<<8)+(((12)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_conn_timing_parameters;
}

/** 
*
* gecko_cmd_le_gap_set_advertise_random_address
*
* Set the advertiser on an advertising set to use a random address. This
* overrides the default advertiser address which is either the public device
* address programmed at production or the address written into persistent
* storage using system_set_identity_address command. This setting is stored in
* RAM only and does not change the identity address in persistent storage.
* 
* When setting a resolvable random address, the address parameter is ignored.
* The stack generates a private resolvable random address and set it as the
* advertiser address. The generated address is returned in the response.
* 
* To use the default advertiser address, remove this setting using
* le_gap_clear_advertise_random_address command.
* 
* Wrong state error is returned if advertising has been enabled on the
* advertising set. Invalid parameter error is returned if the advertising set
* handle is invalid or the address does not conforms to the Bluetooth
* specification. 
*
* @param handle   Advertising set handle
* @param addr_type   Address type:
*  
*      1: Static device address  
*      2: Private resolvable random address  
*      3: Private non-resolvable random address. This type can only be used for
*      non-connectable advertising.
* @param address   The random address to set. Ignore this field when setting a resolvable random
*  address.
*
**/

static inline struct gecko_msg_le_gap_set_advertise_random_address_rsp_t* gecko_cmd_le_gap_set_advertise_random_address(uint8 handle,uint8 addr_type,bd_addr address)
{
    
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_random_address.handle=handle;
    gecko_cmd_msg->data.cmd_le_gap_set_advertise_random_address.addr_type=addr_type;
    memcpy(&gecko_cmd_msg->data.cmd_le_gap_set_advertise_random_address.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->header=(gecko_cmd_le_gap_set_advertise_random_address_id+(((8)&0xff)<<8)+(((8)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_set_advertise_random_address;
}

/** 
*
* gecko_cmd_le_gap_clear_advertise_random_address
*
* Clear the random address previously set for the advertiser address on an
* advertising set. A random address can be set using
* le_gap_set_advertise_random_address command. The default advertiser address
* will be used after this operation.
* 
* Wrong state error is returned if advertising has been enabled on the
* advertising set. Invalid parameter error is returned if the advertising set
* handle is invalid. 
*
* @param handle   Advertising set handle
*
**/

static inline struct gecko_msg_le_gap_clear_advertise_random_address_rsp_t* gecko_cmd_le_gap_clear_advertise_random_address(uint8 handle)
{
    
    gecko_cmd_msg->data.cmd_le_gap_clear_advertise_random_address.handle=handle;
    gecko_cmd_msg->header=(gecko_cmd_le_gap_clear_advertise_random_address_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_gap_clear_advertise_random_address;
}

/** 
*
* gecko_cmd_sync_open
*
* Establish a synchronization with a periodic advertising from the specified
* advertiser and begin receiving periodic advertising packets. Note that
* synchronization establishment can only occur when scanning is enabled. While
* scanning is disabled, no attempt to synchronize will occur.
* 
* The application should determine skip and timeout values based on the periodic
* advertising interval provided by the advertiser. If skip and timeout are used,
* select appropriate values so that they allow a few receiving attempts.
* Periodic advertising intervals are reported in event
* le_gap_extended_scan_response. 
*
* @param adv_sid   Advertising set identifier
* @param skip   The maximum number of periodic advertising packets that can be skipped after a
*  successful receive. Range: 0x0000 to 0x01F3
* @param timeout   The maximum permitted time between successful receives. If this time is
*  exceeded, synchronization is lost. Unit: 10 ms.
*  
*      Range: 0x06 to 0xFFFF  
*      Unit: 10 ms  
*      Time range: 100 ms ms to 163.84 s
* @param address   Address of the advertiser
* @param address_type   Advertiser address type. Values:
*  
*      0: Public address  
*       1: Random address
*
* Events generated
*
* gecko_evt_sync_opened - Triggered after the synchronization is established.
* gecko_evt_sync_data - Indicates that a periodic advertisement packet is received.
*
**/

static inline struct gecko_msg_sync_open_rsp_t* gecko_cmd_sync_open(uint8 adv_sid,uint16 skip,uint16 timeout,bd_addr address,uint8 address_type)
{
    
    gecko_cmd_msg->data.cmd_sync_open.adv_sid=adv_sid;
    gecko_cmd_msg->data.cmd_sync_open.skip=skip;
    gecko_cmd_msg->data.cmd_sync_open.timeout=timeout;
    memcpy(&gecko_cmd_msg->data.cmd_sync_open.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->data.cmd_sync_open.address_type=address_type;
    gecko_cmd_msg->header=(gecko_cmd_sync_open_id+(((12)&0xff)<<8)+(((12)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sync_open;
}

/** 
*
* gecko_cmd_sync_close
*
* Closes a periodic advertising synchronization or cancels an ongoing attempt of
* establishing a synchronization. 
*
* @param sync   Periodic advertising synchronization handle
*
* Events generated
*
* gecko_evt_sync_closed - Triggered after a periodic advertising synchronization has been closed or
*  canceled.
*
**/

static inline struct gecko_msg_sync_close_rsp_t* gecko_cmd_sync_close(uint8 sync)
{
    
    gecko_cmd_msg->data.cmd_sync_close.sync=sync;
    gecko_cmd_msg->header=(gecko_cmd_sync_close_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sync_close;
}

/** 
*
* gecko_cmd_le_connection_set_parameters
*
* Deprecated and replaced by le_connection_set_timing_parameters command.
* 
* Request a change in the connection parameters of a Bluetooth connection. 
*
* @param connection   Connection Handle
* @param min_interval   Minimum value for the connection event interval. This must be set be less than
*  or equal to max_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s
* @param max_interval   Maximum value for the connection event interval. This must be set greater than
*  or equal to min_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s
* @param latency   Slave latency, which defines how many connection intervals the slave can skip
*  if it has no data to send
*  
*      Range: 0x0000 to 0x01f4  
*  
*  Use 0x0000 for default value
* @param timeout   Supervision timeout, which defines the time that the connection is maintained
*  although the devices can't communicate at the currently configured connection
*  intervals.
*  
*      Range: 0x000a to 0x0c80  
*      Time = Value x 10 ms  
*      Time Range: 100 ms to 32 s  
*      The value in milliseconds must be larger than (1 + latency) * max_interval
*      * 2, where max_interval is given in milliseconds  
*  
*  Set the supervision timeout at a value which allows communication attempts
*  over at least a few connection intervals.
*
* Events generated
*
* gecko_evt_le_connection_parameters - Triggered after new connection parameters are applied on the connection.
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_connection_set_parameters_rsp_t* gecko_cmd_le_connection_set_parameters(uint8 connection,uint16 min_interval,uint16 max_interval,uint16 latency,uint16 timeout)
{
    
    gecko_cmd_msg->data.cmd_le_connection_set_parameters.connection=connection;
    gecko_cmd_msg->data.cmd_le_connection_set_parameters.min_interval=min_interval;
    gecko_cmd_msg->data.cmd_le_connection_set_parameters.max_interval=max_interval;
    gecko_cmd_msg->data.cmd_le_connection_set_parameters.latency=latency;
    gecko_cmd_msg->data.cmd_le_connection_set_parameters.timeout=timeout;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_set_parameters_id+(((9)&0xff)<<8)+(((9)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_set_parameters;
}

/** 
*
* gecko_cmd_le_connection_get_rssi
*
* Get the latest RSSI value of a Bluetooth connection. The RSSI value will be
* reported in a le_connection_rssi event. 
*
* @param connection   Connection handle
*
* Events generated
*
* gecko_evt_le_connection_rssi - Triggered when this command has completed.
*
**/

static inline struct gecko_msg_le_connection_get_rssi_rsp_t* gecko_cmd_le_connection_get_rssi(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_le_connection_get_rssi.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_get_rssi_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_get_rssi;
}

/** 
*
* gecko_cmd_le_connection_disable_slave_latency
*
* Temporarily enable or disable slave latency. Used only when Bluetooth device
* is acting as slave. When slave latency is disabled, the slave latency
* connection parameter is not set to 0 but the device will wake up on every
* connection interval to receive and send packets. 
*
* @param connection   Connection Handle
* @param disable   0 enable, 1 disable slave latency. Default: 0
*
**/

static inline struct gecko_msg_le_connection_disable_slave_latency_rsp_t* gecko_cmd_le_connection_disable_slave_latency(uint8 connection,uint8 disable)
{
    
    gecko_cmd_msg->data.cmd_le_connection_disable_slave_latency.connection=connection;
    gecko_cmd_msg->data.cmd_le_connection_disable_slave_latency.disable=disable;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_disable_slave_latency_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_disable_slave_latency;
}

/** 
*
* gecko_cmd_le_connection_set_phy
*
* Deprecated and replaced by le_connection_set_preferred_phy command.
* 
* Set preferred PHYs for a connection. Preferred PHYs are connection-specific.
* Event le_connection_phy_status is received when PHY update procedure is
* completed. Non-preferred PHY can also be set if remote device does not accept
* any of the preferred PHYs.
* 
* NOTE: 2 Mbit and Coded PHYs are not supported by all devices. 
*
* @param connection   
* @param phy   Preferred PHYs for connection. This parameter is a bitfield and multiple PHYs
*  can be preferred by setting multiple bits.
*  
*      0x01: 1M PHY  
*      0x02: 2M PHY  
*      0x04: 125k Coded PHY (S=8)  
*      0x08: 500k Coded PHY (S=2)
*
* Events generated
*
* gecko_evt_le_connection_phy_status - 
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_le_connection_set_phy_rsp_t* gecko_cmd_le_connection_set_phy(uint8 connection,uint8 phy)
{
    
    gecko_cmd_msg->data.cmd_le_connection_set_phy.connection=connection;
    gecko_cmd_msg->data.cmd_le_connection_set_phy.phy=phy;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_set_phy_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_set_phy;
}

/** 
*
* gecko_cmd_le_connection_close
*
* Close a Bluetooth connection or cancel an ongoing connection establishment
* process. The parameter is a connection handle which is reported in
* le_connection_opened event or le_gap_connect response. 
*
* @param connection   Handle of the connection to be closed
*
* Events generated
*
* gecko_evt_le_connection_closed - 
*
**/

static inline struct gecko_msg_le_connection_close_rsp_t* gecko_cmd_le_connection_close(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_le_connection_close.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_close_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_close;
}

/** 
*
* gecko_cmd_le_connection_set_timing_parameters
*
* Request a change in the connection parameters of a Bluetooth connection. 
*
* @param connection   Connection Handle
* @param min_interval   Minimum value for the connection event interval. This must be set less than or
*  equal to max_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s
* @param max_interval   Maximum value for the connection event interval. This must be set greater than
*  or equal to min_interval.
*  
*      Time = Value x 1.25 ms  
*      Range: 0x0006 to 0x0c80  
*      Time Range: 7.5 ms to 4 s
* @param latency   Slave latency, which defines how many connection intervals the slave can skip
*  if it has no data to send
*  
*      Range: 0x0000 to 0x01f4  
*  
*  Use 0x0000 for default value
* @param timeout   Supervision timeout, which defines the time that the connection is maintained
*  although the devices can't communicate at the currently configured connection
*  intervals.
*  
*      Range: 0x000a to 0x0c80  
*      Time = Value x 10 ms  
*      Time Range: 100 ms to 32 s  
*      The value in milliseconds must be larger than (1 + latency) * max_interval
*      * 2, where max_interval is given in milliseconds  
*  
*  Set the supervision timeout at a value which allows communication attempts
*  over at least a few connection intervals.
* @param min_ce_length   Minimum value for the connection event length. This must be set less than or
*  equal to max_ce_length.
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0000 to 0xffff  
*  
*  Value is not currently used and is reserved for future. Set to 0.
* @param max_ce_length   Maximum value for the connection event length. This must be set greater than
*  or equal to min_ce_length.
*  
*      Time = Value x 0.625 ms  
*      Range: 0x0000 to 0xffff  
*  
*  Use 0xffff for no limitation.
*
* Events generated
*
* gecko_evt_le_connection_parameters - Triggered after new connection parameters are applied on the connection.
*
**/

static inline struct gecko_msg_le_connection_set_timing_parameters_rsp_t* gecko_cmd_le_connection_set_timing_parameters(uint8 connection,uint16 min_interval,uint16 max_interval,uint16 latency,uint16 timeout,uint16 min_ce_length,uint16 max_ce_length)
{
    
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.connection=connection;
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.min_interval=min_interval;
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.max_interval=max_interval;
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.latency=latency;
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.timeout=timeout;
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.min_ce_length=min_ce_length;
    gecko_cmd_msg->data.cmd_le_connection_set_timing_parameters.max_ce_length=max_ce_length;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_set_timing_parameters_id+(((13)&0xff)<<8)+(((13)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_set_timing_parameters;
}

/** 
*
* gecko_cmd_le_connection_read_channel_map
*
* Read channel map for a specified connection. 
*
* @param connection   Connection Handle
*
**/

static inline struct gecko_msg_le_connection_read_channel_map_rsp_t* gecko_cmd_le_connection_read_channel_map(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_le_connection_read_channel_map.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_read_channel_map_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_read_channel_map;
}

/** 
*
* gecko_cmd_le_connection_set_preferred_phy
*
* Sets preferred and accepted PHYs for the given connection. Event
* le_connection_phy_status is received when PHY update procedure is completed.
* Non-preferred PHY can also be set if remote device does not accept any of the
* preferred PHYs.
* 
* The parameter accepted_phy is used for specifying the PHYs that the stack can
* accept in a remote initiated PHY update request. A PHY update will not occur
* if none of the accepted PHYs presents in the request.
* 
* NOTE: 2M and Coded PHYs are not supported by all devices. 
*
* @param connection   Connection handle
* @param preferred_phy   Preferred PHYs. This parameter is a bitfield and multiple PHYs can be set.
*  
*      0x01: 1M PHY  
*      0x02: 2M PHY  
*      0x04: 125k Coded PHY (S=8)  
*      0x08: 500k Coded PHY (S=2)  
*  
*  Default: 0xff (no preference)
* @param accepted_phy   Accepted PHYs in remotely-initiated PHY update requests. This parameter is a
*  bitfield and multiple PHYs can be set.
*  
*      0x01: 1M PHY  
*      0x02: 2M PHY  
*      0x04: Coded PHY  
*      0xff: Any PHYs  
*  
*  Default: 0xff (all PHYs accepted)
*
* Events generated
*
* gecko_evt_le_connection_phy_status - 
*
**/

static inline struct gecko_msg_le_connection_set_preferred_phy_rsp_t* gecko_cmd_le_connection_set_preferred_phy(uint8 connection,uint8 preferred_phy,uint8 accepted_phy)
{
    
    gecko_cmd_msg->data.cmd_le_connection_set_preferred_phy.connection=connection;
    gecko_cmd_msg->data.cmd_le_connection_set_preferred_phy.preferred_phy=preferred_phy;
    gecko_cmd_msg->data.cmd_le_connection_set_preferred_phy.accepted_phy=accepted_phy;
    gecko_cmd_msg->header=(gecko_cmd_le_connection_set_preferred_phy_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_le_connection_set_preferred_phy;
}

/** 
*
* gecko_cmd_gatt_set_max_mtu
*
* Set the maximum size of ATT Message Transfer Units (MTU). Functionality is the
* same as gatt_server_set_max_mtu and this setting applies to both GATT client
* and server. If the given value is too large according to the maximum BGAPI
* payload size, the system will select the maximum possible value as the maximum
* ATT_MTU. If maximum ATT_MTU is larger than 23, the GATT client in the stack
* will automatically send an MTU exchange request after a Bluetooth connection
* has been established. 
*
* @param max_mtu   Maximum size of Message Transfer Units (MTU) allowed
*  
*      Range: 23 to 250  
*  
*  Default: 247
*
**/

static inline struct gecko_msg_gatt_set_max_mtu_rsp_t* gecko_cmd_gatt_set_max_mtu(uint16 max_mtu)
{
    
    gecko_cmd_msg->data.cmd_gatt_set_max_mtu.max_mtu=max_mtu;
    gecko_cmd_msg->header=(gecko_cmd_gatt_set_max_mtu_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_set_max_mtu;
}

/** 
*
* gecko_cmd_gatt_discover_primary_services
*
* Discover all primary services of a remote GATT database. This command
* generates a unique gatt_service event for every discovered primary service.
* Received gatt_procedure_completed event indicates that this GATT procedure has
* successfully completed or failed with an error. 
*
* @param connection   
*
* Events generated
*
* gecko_evt_gatt_service - Discovered service from remote GATT database
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_discover_primary_services_rsp_t* gecko_cmd_gatt_discover_primary_services(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_gatt_discover_primary_services.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_gatt_discover_primary_services_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_discover_primary_services;
}

/** 
*
* gecko_cmd_gatt_discover_primary_services_by_uuid
*
* Discover primary services with the specified UUID in a remote GATT database.
* This command generates unique gatt_service event for every discovered primary
* service. Received gatt_procedure_completed event indicates that this GATT
* procedure was successfully completed or failed with an error. 
*
* @param connection   
* @param uuid_len   Array length
* @param uuid_data   Service UUID in little endian format
*
* Events generated
*
* gecko_evt_gatt_service - Discovered service from remote GATT database.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_discover_primary_services_by_uuid_rsp_t* gecko_cmd_gatt_discover_primary_services_by_uuid(uint8 connection,uint8 uuid_len, const uint8* uuid_data)
{
    if ((uint16_t)uuid_len > BGLIB_MSG_MAX_PAYLOAD - 2)
    {
        gecko_rsp_msg->data.rsp_gatt_discover_primary_services_by_uuid.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_discover_primary_services_by_uuid;
    }

    
    gecko_cmd_msg->data.cmd_gatt_discover_primary_services_by_uuid.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_discover_primary_services_by_uuid.uuid.len=uuid_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_discover_primary_services_by_uuid.uuid.data,uuid_data,uuid_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_discover_primary_services_by_uuid_id+(((2+uuid_len)&0xff)<<8)+(((2+uuid_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_discover_primary_services_by_uuid;
}

/** 
*
* gecko_cmd_gatt_discover_characteristics
*
* Discover all characteristics of a GATT service from a remote GATT database.
* This command generates a unique gatt_characteristic event for every discovered
* characteristic. Received gatt_procedure_completed event indicates that this
* GATT procedure was successfully completed or failed with an error. 
*
* @param connection   
* @param service   
*
* Events generated
*
* gecko_evt_gatt_characteristic - Discovered characteristic from remote GATT database.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_discover_characteristics_rsp_t* gecko_cmd_gatt_discover_characteristics(uint8 connection,uint32 service)
{
    
    gecko_cmd_msg->data.cmd_gatt_discover_characteristics.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_discover_characteristics.service=service;
    gecko_cmd_msg->header=(gecko_cmd_gatt_discover_characteristics_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_discover_characteristics;
}

/** 
*
* gecko_cmd_gatt_discover_characteristics_by_uuid
*
* Discover all characteristics of a GATT service in a remote GATT database
* having the specified UUID. This command generates a unique gatt_characteristic
* event for every discovered characteristic having the specified UUID. Received
* gatt_procedure_completed event indicates that this GATT procedure was
* successfully completed or failed with an error. 
*
* @param connection   
* @param service   
* @param uuid_len   Array length
* @param uuid_data   Characteristic UUID in little endian format
*
* Events generated
*
* gecko_evt_gatt_characteristic - Discovered characteristic from remote GATT database.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_discover_characteristics_by_uuid_rsp_t* gecko_cmd_gatt_discover_characteristics_by_uuid(uint8 connection,uint32 service,uint8 uuid_len, const uint8* uuid_data)
{
    if ((uint16_t)uuid_len > BGLIB_MSG_MAX_PAYLOAD - 6)
    {
        gecko_rsp_msg->data.rsp_gatt_discover_characteristics_by_uuid.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_discover_characteristics_by_uuid;
    }

    
    gecko_cmd_msg->data.cmd_gatt_discover_characteristics_by_uuid.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_discover_characteristics_by_uuid.service=service;
    gecko_cmd_msg->data.cmd_gatt_discover_characteristics_by_uuid.uuid.len=uuid_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_discover_characteristics_by_uuid.uuid.data,uuid_data,uuid_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_discover_characteristics_by_uuid_id+(((6+uuid_len)&0xff)<<8)+(((6+uuid_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_discover_characteristics_by_uuid;
}

/** 
*
* gecko_cmd_gatt_set_characteristic_notification
*
* Enable or disable the notifications and indications sent from a remote GATT
* server. This procedure discovers a characteristic client configuration
* descriptor and writes the related configuration flags to a remote GATT
* database. A received gatt_procedure_completed event indicates that this GATT
* procedure was successfully completed or that it failed with an error. 
*
* @param connection   
* @param characteristic   
* @param flags   Characteristic client configuration flags
*
* Events generated
*
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
* gecko_evt_gatt_characteristic_value - If an indication or notification has been enabled for a characteristic, this
*  event is triggered whenever an indication or notification is sent by the
*  remote GATT server. The triggering conditions of the GATT server are defined
*  by an upper level, for example by a profile. As a result, it is possible that
*  no values are ever received, or that it may take time, depending on how the
*  server is configured.
*
**/

static inline struct gecko_msg_gatt_set_characteristic_notification_rsp_t* gecko_cmd_gatt_set_characteristic_notification(uint8 connection,uint16 characteristic,uint8 flags)
{
    
    gecko_cmd_msg->data.cmd_gatt_set_characteristic_notification.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_set_characteristic_notification.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_set_characteristic_notification.flags=flags;
    gecko_cmd_msg->header=(gecko_cmd_gatt_set_characteristic_notification_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_set_characteristic_notification;
}

/** 
*
* gecko_cmd_gatt_discover_descriptors
*
* Discover all descriptors of a GATT characteristic in a remote GATT database.
* It generates a unique gatt_descriptor event for every discovered descriptor.
* Received gatt_procedure_completed event indicates that this GATT procedure has
* successfully completed or failed with an error. 
*
* @param connection   
* @param characteristic   
*
* Events generated
*
* gecko_evt_gatt_descriptor - Discovered descriptor from remote GATT database.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_discover_descriptors_rsp_t* gecko_cmd_gatt_discover_descriptors(uint8 connection,uint16 characteristic)
{
    
    gecko_cmd_msg->data.cmd_gatt_discover_descriptors.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_discover_descriptors.characteristic=characteristic;
    gecko_cmd_msg->header=(gecko_cmd_gatt_discover_descriptors_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_discover_descriptors;
}

/** 
*
* gecko_cmd_gatt_read_characteristic_value
*
* Read the value of a characteristic from a remote GATT database. A single
* gatt_characteristic_value event is generated if the characteristic value fits
* in one ATT PDU. Otherwise, more than one  gatt_characteristic_value event is
* generated because the firmware will automatically use the Read Long
* Characteristic Values procedure. A received gatt_procedure_completed event
* indicates that all data was read successfully or that an error response was
* received.
* 
* Note that the GATT client does not verify if the requested attribute is a
* characteristic value. Therefore, before calling this command, ensure that the
* attribute handle is for a characteristic value, for example, by performing
* characteristic discovery. 
*
* @param connection   
* @param characteristic   
*
* Events generated
*
* gecko_evt_gatt_characteristic_value - Contains the data of a characteristic sent by the GATT Server.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_read_characteristic_value_rsp_t* gecko_cmd_gatt_read_characteristic_value(uint8 connection,uint16 characteristic)
{
    
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value.characteristic=characteristic;
    gecko_cmd_msg->header=(gecko_cmd_gatt_read_characteristic_value_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_read_characteristic_value;
}

/** 
*
* gecko_cmd_gatt_read_characteristic_value_by_uuid
*
* Read characteristic values of a service from a remote GATT database by giving
* the UUID of the characteristic and the handle of the service containing this
* characteristic. If multiple characteristic values are received in one ATT PDU,
* then one  gatt_characteristic_value event is generated for each value. If the
* first characteristic value does not fit in one ATT PDU, the firmware
* automatically uses the Read Long Characteristic Values procedure and generate
* more  gatt_characteristic_value events until the value has been completely
* read. A received gatt_procedure_completed event indicates that all data was
* read successfully or that an error response was received. 
*
* @param connection   
* @param service   
* @param uuid_len   Array length
* @param uuid_data   Characteristic UUID in little endian format
*
* Events generated
*
* gecko_evt_gatt_characteristic_value - Contains the data of a characteristic sent by the GATT Server.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_read_characteristic_value_by_uuid_rsp_t* gecko_cmd_gatt_read_characteristic_value_by_uuid(uint8 connection,uint32 service,uint8 uuid_len, const uint8* uuid_data)
{
    if ((uint16_t)uuid_len > BGLIB_MSG_MAX_PAYLOAD - 6)
    {
        gecko_rsp_msg->data.rsp_gatt_read_characteristic_value_by_uuid.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_read_characteristic_value_by_uuid;
    }

    
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_by_uuid.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_by_uuid.service=service;
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_by_uuid.uuid.len=uuid_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_by_uuid.uuid.data,uuid_data,uuid_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_read_characteristic_value_by_uuid_id+(((6+uuid_len)&0xff)<<8)+(((6+uuid_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_read_characteristic_value_by_uuid;
}

/** 
*
* gecko_cmd_gatt_write_characteristic_value
*
* Write the value of a characteristic in a remote GATT database. If the given
* value does not fit in one ATT PDU, "write long" GATT procedure is used
* automatically. Received gatt_procedure_completed event indicates that all data
* was written successfully or that an error response was received. 
*
* @param connection   
* @param characteristic   
* @param value_len   Array length
* @param value_data   Characteristic value
*
* Events generated
*
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_write_characteristic_value_rsp_t* gecko_cmd_gatt_write_characteristic_value(uint8 connection,uint16 characteristic,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_gatt_write_characteristic_value.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_write_characteristic_value;
    }

    
    gecko_cmd_msg->data.cmd_gatt_write_characteristic_value.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_write_characteristic_value.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_write_characteristic_value.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_write_characteristic_value.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_write_characteristic_value_id+(((4+value_len)&0xff)<<8)+(((4+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_write_characteristic_value;
}

/** 
*
* gecko_cmd_gatt_write_characteristic_value_without_response
*
* Write the value of a characteristic in a remote GATT server. It does not
* generate an event. All failures on the server are ignored silently. For
* example, if an error is generated in the remote GATT server and the given
* value is not written into the database, no error message will be reported to
* the local GATT client. Note that this command can't be used to write long
* values. At most ATT_MTU - 3 amount of data can be sent once. 
*
* @param connection   
* @param characteristic   
* @param value_len   Array length
* @param value_data   Characteristic value
*
**/

static inline struct gecko_msg_gatt_write_characteristic_value_without_response_rsp_t* gecko_cmd_gatt_write_characteristic_value_without_response(uint8 connection,uint16 characteristic,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_gatt_write_characteristic_value_without_response.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_write_characteristic_value_without_response;
    }

    
    gecko_cmd_msg->data.cmd_gatt_write_characteristic_value_without_response.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_write_characteristic_value_without_response.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_write_characteristic_value_without_response.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_write_characteristic_value_without_response.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_write_characteristic_value_without_response_id+(((4+value_len)&0xff)<<8)+(((4+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_write_characteristic_value_without_response;
}

/** 
*
* gecko_cmd_gatt_prepare_characteristic_value_write
*
* Add a characteristic value to the write queue of a remote GATT server. It can
* be used when long attributes need to be written or a set of values needs to be
* written atomically. At most ATT_MTU - 5 amount of data can be sent at one
* time. Writes are executed or canceled with the
* execute_characteristic_value_write command. Whether the writes succeed or not
* is indicated in the response of the execute_characteristic_value_write
* command.
* 
* In all use cases where the amount of data to transfer fits into the BGAPI
* payload, use the command gatt_write_characteristic_value to write long values
* because it transparently performs the prepare_write and execute_write
* commands. 
*
* @param connection   
* @param characteristic   
* @param offset   Offset of the characteristic value
* @param value_len   Array length
* @param value_data   Value to write into the specified characteristic of the remote GATT database
*
* Events generated
*
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_prepare_characteristic_value_write_rsp_t* gecko_cmd_gatt_prepare_characteristic_value_write(uint8 connection,uint16 characteristic,uint16 offset,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 6)
    {
        gecko_rsp_msg->data.rsp_gatt_prepare_characteristic_value_write.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_prepare_characteristic_value_write;
    }

    
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_write.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_write.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_write.offset=offset;
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_write.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_write.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_prepare_characteristic_value_write_id+(((6+value_len)&0xff)<<8)+(((6+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_prepare_characteristic_value_write;
}

/** 
*
* gecko_cmd_gatt_execute_characteristic_value_write
*
* Commit or cancel previously queued writes to a long characteristic of a remote
* GATT server. Writes are sent to the queue with
* prepare_characteristic_value_write command. Content, offset, and length of
* queued values are validated by this procedure. A received
* gatt_procedure_completed event indicates that all data was written
* successfully or that an error response was received. 
*
* @param connection   
* @param flags   gatt_execute_write_flag
*  
*      0: cancel  
*      1: commit
*
* Events generated
*
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_execute_characteristic_value_write_rsp_t* gecko_cmd_gatt_execute_characteristic_value_write(uint8 connection,uint8 flags)
{
    
    gecko_cmd_msg->data.cmd_gatt_execute_characteristic_value_write.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_execute_characteristic_value_write.flags=flags;
    gecko_cmd_msg->header=(gecko_cmd_gatt_execute_characteristic_value_write_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_execute_characteristic_value_write;
}

/** 
*
* gecko_cmd_gatt_send_characteristic_confirmation
*
* Send a confirmation to a remote GATT server after receiving a characteristic
* indication. The gatt_characteristic_value event carries the att_opcode
* containing handle_value_indication (0x1d), which reveals that an indication
* has been received and must be confirmed with this command. The confirmation
* needs to be sent within 30 seconds, otherwise further GATT transactions are
* not allowed by the remote side. 
*
* @param connection   
*
**/

static inline struct gecko_msg_gatt_send_characteristic_confirmation_rsp_t* gecko_cmd_gatt_send_characteristic_confirmation(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_gatt_send_characteristic_confirmation.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_gatt_send_characteristic_confirmation_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_send_characteristic_confirmation;
}

/** 
*
* gecko_cmd_gatt_read_descriptor_value
*
* Read the descriptor value of a characteristic in a remote GATT database. A
* single  gatt_descriptor_value event is generated if the descriptor value fits
* in one ATT PDU. Otherwise, more than one gatt_descriptor_value events are
* generated because the firmware automatically uses the Read Long Characteristic
* Values procedure. A received gatt_procedure_completed event indicates that all
* data was read successfully or that an error response was received. 
*
* @param connection   
* @param descriptor   
*
* Events generated
*
* gecko_evt_gatt_descriptor_value - Descriptor value received from the remote GATT server.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_read_descriptor_value_rsp_t* gecko_cmd_gatt_read_descriptor_value(uint8 connection,uint16 descriptor)
{
    
    gecko_cmd_msg->data.cmd_gatt_read_descriptor_value.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_read_descriptor_value.descriptor=descriptor;
    gecko_cmd_msg->header=(gecko_cmd_gatt_read_descriptor_value_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_read_descriptor_value;
}

/** 
*
* gecko_cmd_gatt_write_descriptor_value
*
* Write the value of a characteristic descriptor in a remote GATT database. If
* the given value does not fit in one ATT PDU, "write long" GATT procedure is
* used automatically. Received gatt_procedure_completed event indicates either
* that all data was written successfully or that an error response was received. 
*
* @param connection   
* @param descriptor   
* @param value_len   Array length
* @param value_data   Descriptor value
*
* Events generated
*
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_write_descriptor_value_rsp_t* gecko_cmd_gatt_write_descriptor_value(uint8 connection,uint16 descriptor,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_gatt_write_descriptor_value.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_write_descriptor_value;
    }

    
    gecko_cmd_msg->data.cmd_gatt_write_descriptor_value.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_write_descriptor_value.descriptor=descriptor;
    gecko_cmd_msg->data.cmd_gatt_write_descriptor_value.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_write_descriptor_value.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_write_descriptor_value_id+(((4+value_len)&0xff)<<8)+(((4+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_write_descriptor_value;
}

/** 
*
* gecko_cmd_gatt_find_included_services
*
* Find the services that are included by a service in a remote GATT database.
* This command generates a unique gatt_service event for each included service.
* The received gatt_procedure_completed event indicates that this GATT procedure
* was successfully completed or failed with an error. 
*
* @param connection   
* @param service   
*
* Events generated
*
* gecko_evt_gatt_service - Discovered service from remote GATT database.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_find_included_services_rsp_t* gecko_cmd_gatt_find_included_services(uint8 connection,uint32 service)
{
    
    gecko_cmd_msg->data.cmd_gatt_find_included_services.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_find_included_services.service=service;
    gecko_cmd_msg->header=(gecko_cmd_gatt_find_included_services_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_find_included_services;
}

/** 
*
* gecko_cmd_gatt_read_multiple_characteristic_values
*
* Read values of multiple characteristics from a remote GATT database at once.
* The GATT server returns values in one ATT PDU as the response. If the total
* set of values is greater than (ATT_MTU - 1) bytes in length, only the first
* (ATT_MTU - 1) bytes are included. A single gatt_characteristic_value event is
* generated, in which the characteristic is set to 0 and data in the value
* parameter is a concatenation of characteristic values in the order they were
* requested. The received gatt_procedure_completed event indicates either that
* this GATT procedure was successfully completed or failed with an error.
* 
* Use this command only for characteristics values that have a known fixed size,
* except the last one that could have variable length.
* 
* When the remote GATT server is from Silicon Labs Bluetooth stack, the server
* returns ATT Invalid PDU (0x04) if this command only reads one characteristic
* value. The server returns ATT Application Error 0x80 if this command reads the
* value of a user-type characteristic. 
*
* @param connection   
* @param characteristic_list_len   Array length
* @param characteristic_list_data   List of uint16 characteristic handles each in little endian format.
*
* Events generated
*
* gecko_evt_gatt_characteristic_value - A concatenation of characteristic values in the order they were requested
* gecko_evt_gatt_procedure_completed - Procedure was either successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_read_multiple_characteristic_values_rsp_t* gecko_cmd_gatt_read_multiple_characteristic_values(uint8 connection,uint8 characteristic_list_len, const uint8* characteristic_list_data)
{
    if ((uint16_t)characteristic_list_len > BGLIB_MSG_MAX_PAYLOAD - 2)
    {
        gecko_rsp_msg->data.rsp_gatt_read_multiple_characteristic_values.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_read_multiple_characteristic_values;
    }

    
    gecko_cmd_msg->data.cmd_gatt_read_multiple_characteristic_values.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_read_multiple_characteristic_values.characteristic_list.len=characteristic_list_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_read_multiple_characteristic_values.characteristic_list.data,characteristic_list_data,characteristic_list_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_read_multiple_characteristic_values_id+(((2+characteristic_list_len)&0xff)<<8)+(((2+characteristic_list_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_read_multiple_characteristic_values;
}

/** 
*
* gecko_cmd_gatt_read_characteristic_value_from_offset
*
* Read a partial characteristic value with a specified offset and maximum length
* from a remote GATT database. It is equivalent to
* gatt_read_characteristic_value if both the offset and maximum length
* parameters are 0. A single gatt_characteristic_value event is generated if the
* value to read fits in one ATT PDU. Otherwise, more than one
* gatt_characteristic_value events are generated because the firmware will
* automatically use the Read Long Characteristic Values procedure. A received
* gatt_procedure_completed event indicates that all data was read successfully
* or that an error response was received. 
*
* @param connection   
* @param characteristic   
* @param offset   Offset of the characteristic value
* @param maxlen   Maximum bytes to read. If this parameter is 0, all characteristic values
*  starting at a given offset will be read.
*
* Events generated
*
* gecko_evt_gatt_characteristic_value - Contains the data of a characteristic sent by the GATT Server.
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_read_characteristic_value_from_offset_rsp_t* gecko_cmd_gatt_read_characteristic_value_from_offset(uint8 connection,uint16 characteristic,uint16 offset,uint16 maxlen)
{
    
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_from_offset.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_from_offset.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_from_offset.offset=offset;
    gecko_cmd_msg->data.cmd_gatt_read_characteristic_value_from_offset.maxlen=maxlen;
    gecko_cmd_msg->header=(gecko_cmd_gatt_read_characteristic_value_from_offset_id+(((7)&0xff)<<8)+(((7)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_read_characteristic_value_from_offset;
}

/** 
*
* gecko_cmd_gatt_prepare_characteristic_value_reliable_write
*
* Add a characteristic value to the write queue of a remote GATT server and
* verifies whether the value was correctly received by the server. Received
* gatt_procedure_completed event indicates that this GATT procedure was
* successfully completed or failed with an error. Specifically, error code
* 0x0194 (data_corrupted) will be returned if the value received from the GATT
* server's response fails to pass the reliable write verification. At most
* ATT_MTU - 5 amount of data can be sent at one time. Writes are executed or
* canceled with the execute_characteristic_value_write command. Whether the
* writes succeed or not is indicated in the response of the
* execute_characteristic_value_write command. 
*
* @param connection   
* @param characteristic   
* @param offset   Offset of the characteristic value
* @param value_len   Array length
* @param value_data   Value to write into the specified characteristic of the remote GATT database
*
* Events generated
*
* gecko_evt_gatt_procedure_completed - Procedure was successfully completed or failed with an error.
*
**/

static inline struct gecko_msg_gatt_prepare_characteristic_value_reliable_write_rsp_t* gecko_cmd_gatt_prepare_characteristic_value_reliable_write(uint8 connection,uint16 characteristic,uint16 offset,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 6)
    {
        gecko_rsp_msg->data.rsp_gatt_prepare_characteristic_value_reliable_write.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_prepare_characteristic_value_reliable_write;
    }

    
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_reliable_write.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_reliable_write.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_reliable_write.offset=offset;
    gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_reliable_write.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_prepare_characteristic_value_reliable_write.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_prepare_characteristic_value_reliable_write_id+(((6+value_len)&0xff)<<8)+(((6+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_prepare_characteristic_value_reliable_write;
}

/** 
*
* gecko_cmd_gatt_server_read_attribute_value
*
* Read the value of an attribute from a local GATT database. Only (maximum BGAPI
* payload size - 3) amount of data can be read at once. The application can
* continue reading with increased offset value if it receives (maximum BGAPI
* payload size - 3) amount of data. 
*
* @param attribute   Attribute handle
* @param offset   Value offset
*
**/

static inline struct gecko_msg_gatt_server_read_attribute_value_rsp_t* gecko_cmd_gatt_server_read_attribute_value(uint16 attribute,uint16 offset)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_read_attribute_value.attribute=attribute;
    gecko_cmd_msg->data.cmd_gatt_server_read_attribute_value.offset=offset;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_read_attribute_value_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_read_attribute_value;
}

/** 
*
* gecko_cmd_gatt_server_read_attribute_type
*
* Read the type of an attribute from a local GATT database. The type is a UUID,
* usually 16 or 128 bits long in little endian format. 
*
* @param attribute   Attribute handle
*
**/

static inline struct gecko_msg_gatt_server_read_attribute_type_rsp_t* gecko_cmd_gatt_server_read_attribute_type(uint16 attribute)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_read_attribute_type.attribute=attribute;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_read_attribute_type_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_read_attribute_type;
}

/** 
*
* gecko_cmd_gatt_server_write_attribute_value
*
* Write the value of an attribute in the local GATT database. Writing the value
* of a characteristic of the local GATT database will not trigger notifications
* or indications to the remote GATT client if the characteristic has a property
* to indicate or notify and the client has enabled notification or indication.
* Notifications and indications are sent to the remote GATT client using
* gatt_server_send_characteristic_notification command. 
*
* @param attribute   Attribute handle
* @param offset   Value offset
* @param value_len   Array length
* @param value_data   Value
*
**/

static inline struct gecko_msg_gatt_server_write_attribute_value_rsp_t* gecko_cmd_gatt_server_write_attribute_value(uint16 attribute,uint16 offset,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 5)
    {
        gecko_rsp_msg->data.rsp_gatt_server_write_attribute_value.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_server_write_attribute_value;
    }

    
    gecko_cmd_msg->data.cmd_gatt_server_write_attribute_value.attribute=attribute;
    gecko_cmd_msg->data.cmd_gatt_server_write_attribute_value.offset=offset;
    gecko_cmd_msg->data.cmd_gatt_server_write_attribute_value.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_server_write_attribute_value.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_write_attribute_value_id+(((5+value_len)&0xff)<<8)+(((5+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_write_attribute_value;
}

/** 
*
* gecko_cmd_gatt_server_send_user_read_response
*
* Send a response to a user_read_request event. The response needs to be sent
* within 30 seconds, otherwise no more GATT transactions are allowed by the
* remote side. If attr_errorcode is set to 0, the characteristic value is sent
* to the remote GATT client in the standard way. Other attr_errorcode values
* will cause the local GATT server to send an attribute protocol error response
* instead of the actual data. At most ATT_MTU - 1 amount of data can be sent at
* one time. The client will continue reading by sending new read request with an
* increased offset value if it receives ATT_MTU - 1 amount of data. 
*
* @param connection   
* @param characteristic   
* @param att_errorcode   
* @param value_len   Array length
* @param value_data   Characteristic value to send to the GATT client. Ignored if att_errorcode is
*  not 0.
*
**/

static inline struct gecko_msg_gatt_server_send_user_read_response_rsp_t* gecko_cmd_gatt_server_send_user_read_response(uint8 connection,uint16 characteristic,uint8 att_errorcode,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 5)
    {
        gecko_rsp_msg->data.rsp_gatt_server_send_user_read_response.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_server_send_user_read_response;
    }

    
    gecko_cmd_msg->data.cmd_gatt_server_send_user_read_response.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_server_send_user_read_response.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_server_send_user_read_response.att_errorcode=att_errorcode;
    gecko_cmd_msg->data.cmd_gatt_server_send_user_read_response.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_server_send_user_read_response.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_send_user_read_response_id+(((5+value_len)&0xff)<<8)+(((5+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_send_user_read_response;
}

/** 
*
* gecko_cmd_gatt_server_send_user_write_response
*
* Send a response to a gatt_server_user_write_request event when parameter
* att_opcode in the event is Write Request (see att_opcode). The response needs
* to be sent within 30 seconds, otherwise no more GATT transactions are allowed
* by the remote side. If attr_errorcode is set to 0, the ATT protocol's write
* response is sent to indicate to the remote GATT client that the write
* operation was processed successfully. Other values will cause the local GATT
* server to send an ATT protocol error response. 
*
* @param connection   
* @param characteristic   
* @param att_errorcode   
*
**/

static inline struct gecko_msg_gatt_server_send_user_write_response_rsp_t* gecko_cmd_gatt_server_send_user_write_response(uint8 connection,uint16 characteristic,uint8 att_errorcode)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_send_user_write_response.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_server_send_user_write_response.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_server_send_user_write_response.att_errorcode=att_errorcode;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_send_user_write_response_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_send_user_write_response;
}

/** 
*
* gecko_cmd_gatt_server_send_characteristic_notification
*
* Send notifications or indications to one or more remote GATT clients. At most
* ATT_MTU - 3 amount of data can be sent one time.
* 
* A notification or indication is sent only if the client has enabled it by
* setting the corresponding flag to the Client Characteristic Configuration
* descriptor. If the Client Characteristic Configuration descriptor supports
* both notifications and indications, the stack will always send a notification
* even when the client has enabled both.
* 
* A new indication to a GATT client can't be sent until an outstanding
* indication procedure with the same client has completed. The procedure is
* completed when a confirmation from the client is received. The confirmation is
* indicated by gatt_server_characteristic_status event.
* 
* Error bg_err_wrong_state is returned if the characteristic does not have the
* notification property, or if the client has not enabled the notification. The
* same applies to the indication property, and in addition, bg_err_wrong_state
* is returned if an indication procedure with the same client is outstanding.
* Always check the response for this command for errors before trying to send
* more data. 
*
* @param connection   A handle of the connection over which the notification or indication is sent.
*  Values:
*  
*      0xff: Sends notification or indication to all connected devices.  
*       Other: Connection handle
* @param characteristic   Characteristic handle
* @param value_len   Array length
* @param value_data   Value to be notified or indicated
*
**/

static inline struct gecko_msg_gatt_server_send_characteristic_notification_rsp_t* gecko_cmd_gatt_server_send_characteristic_notification(uint8 connection,uint16 characteristic,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_gatt_server_send_characteristic_notification.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_server_send_characteristic_notification;
    }

    
    gecko_cmd_msg->data.cmd_gatt_server_send_characteristic_notification.connection=connection;
    gecko_cmd_msg->data.cmd_gatt_server_send_characteristic_notification.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_gatt_server_send_characteristic_notification.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_server_send_characteristic_notification.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_send_characteristic_notification_id+(((4+value_len)&0xff)<<8)+(((4+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_send_characteristic_notification;
}

/** 
*
* gecko_cmd_gatt_server_find_attribute
*
* Find attributes of a certain type from a local GATT database. The type is
* usually given as a 16-bit or 128-bit UUID in little endian format. 
*
* @param start   Search start handle
* @param type_len   Array length
* @param type_data   
*
**/

static inline struct gecko_msg_gatt_server_find_attribute_rsp_t* gecko_cmd_gatt_server_find_attribute(uint16 start,uint8 type_len, const uint8* type_data)
{
    if ((uint16_t)type_len > BGLIB_MSG_MAX_PAYLOAD - 3)
    {
        gecko_rsp_msg->data.rsp_gatt_server_find_attribute.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_gatt_server_find_attribute;
    }

    
    gecko_cmd_msg->data.cmd_gatt_server_find_attribute.start=start;
    gecko_cmd_msg->data.cmd_gatt_server_find_attribute.type.len=type_len;
    memcpy(gecko_cmd_msg->data.cmd_gatt_server_find_attribute.type.data,type_data,type_len);
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_find_attribute_id+(((3+type_len)&0xff)<<8)+(((3+type_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_find_attribute;
}

/** 
*
* gecko_cmd_gatt_server_set_capabilities
*
* Reset capabilities that should be enabled by the GATT database. A service is
* visible to remote GATT clients if at least one of its capabilities is enabled.
* The same applies to a characteristic and its attributes. Capability
* identifiers and their corresponding bit flag values can be found in the auto-
* generated database header file. See UG118: Blue Gecko Bluetooth Profile
* Toolkit Developer's Guide for how to declare capabilities in the GATT
* database.
* 
* Changing the capabilities of a database effectively causes a database change
* (attributes being added or removed) from a remote GATT client point of view.
* If the database has a Generic Attribute service and Service Changed
* characteristic, the stack will monitor the local database change status and
* manage service changed indications for a GATT client that has enabled the
* indication configuration of the Service Changed characteristic. 
*
* @param caps   Bit flags of capabilities to reset. Value 0 sets the default database
*  capabilities.
* @param reserved   Use the value 0 on this reserved field. Do not use none-zero values because
*  they are reserved for future use.
*
**/

static inline struct gecko_msg_gatt_server_set_capabilities_rsp_t* gecko_cmd_gatt_server_set_capabilities(uint32 caps,uint32 reserved)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_set_capabilities.caps=caps;
    gecko_cmd_msg->data.cmd_gatt_server_set_capabilities.reserved=reserved;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_set_capabilities_id+(((8)&0xff)<<8)+(((8)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_set_capabilities;
}

/** 
*
* gecko_cmd_gatt_server_set_max_mtu
*
* Set the maximum size of ATT Message Transfer Units (MTU). The functionality is
* the same as gatt_set_max_mtu and this setting applies to both GATT client and
* server. If the given value is too large according to the maximum BGAPI payload
* size, the system will select the maximum possible value as the maximum
* ATT_MTU. If the maximum ATT_MTU is larger than 23, the GATT client in the
* stack will automatically send an MTU exchange request after a Bluetooth
* connection was established. 
*
* @param max_mtu   Maximum size of Message Transfer Units (MTU) allowed
*  
*      Range: 23 to 250  
*  
*  Default: 247
*
**/

static inline struct gecko_msg_gatt_server_set_max_mtu_rsp_t* gecko_cmd_gatt_server_set_max_mtu(uint16 max_mtu)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_set_max_mtu.max_mtu=max_mtu;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_set_max_mtu_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_set_max_mtu;
}

/** 
*
* gecko_cmd_gatt_server_get_mtu
*
* Get the size of ATT Message Transfer Units (MTU) for a connection. 
*
* @param connection   
*
**/

static inline struct gecko_msg_gatt_server_get_mtu_rsp_t* gecko_cmd_gatt_server_get_mtu(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_get_mtu.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_get_mtu_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_get_mtu;
}

/** 
*
* gecko_cmd_gatt_server_enable_capabilities
*
* Enable additional capabilities in the local GATT database. Already enabled
* capabilities keep unchanged after this command. See
* gatt_server_set_capabilities for more formation. 
*
* @param caps   Capabilities to enable
*
**/

static inline struct gecko_msg_gatt_server_enable_capabilities_rsp_t* gecko_cmd_gatt_server_enable_capabilities(uint32 caps)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_enable_capabilities.caps=caps;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_enable_capabilities_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_enable_capabilities;
}

/** 
*
* gecko_cmd_gatt_server_disable_capabilities
*
* Disable the given capabilities in the local GATT database. See
* gatt_server_set_capabilities for more formation. 
*
* @param caps   Capabilities to disable
*
**/

static inline struct gecko_msg_gatt_server_disable_capabilities_rsp_t* gecko_cmd_gatt_server_disable_capabilities(uint32 caps)
{
    
    gecko_cmd_msg->data.cmd_gatt_server_disable_capabilities.caps=caps;
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_disable_capabilities_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_disable_capabilities;
}

/** 
*
* gecko_cmd_gatt_server_get_enabled_capabilities
*
* Get capabilities currently enabled in the local GATT database. 
*
*
**/

static inline struct gecko_msg_gatt_server_get_enabled_capabilities_rsp_t* gecko_cmd_gatt_server_get_enabled_capabilities()
{
    
    gecko_cmd_msg->header=(gecko_cmd_gatt_server_get_enabled_capabilities_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_gatt_server_get_enabled_capabilities;
}

/** 
*
* gecko_cmd_hardware_set_soft_timer
*
* Start a software timer. Multiple concurrent timers can be running
* simultaneously. 256 unique timer handles (IDs) are available. The maximum
* number of concurrent timers is configurable at device initialization. Up to 16
* concurrent timers can be configured. The default configuration is 4. As the
* RAM for storing timer data is pre-allocated at initialization, an application
* should not configure the amount more than it needs for minimizing RAM usage. 
*
* @param time   Frequency interval of events, which indicates how often to send events in
*  hardware clock ticks (1 second is equal to 32768 ticks).
*  
*  The smallest interval value supported is 328, which is around 10 milliseconds.
*  Any parameters between 0 and 328 will be rounded up to 328. The maximum value
*  is 2147483647, which corresponds to about 18.2 hours.
*  
*  If time is 0, removes the scheduled timer with the same handle.
* @param handle   Timer handle to use, which is returned in timeout event
* @param single_shot   Timer mode. Values:
*  
*      0: false (timer is repeating)  
*       1: true (timer runs only once)
*
* Events generated
*
* gecko_evt_hardware_soft_timer - Sent after this timer has lapsed.
*
**/

static inline struct gecko_msg_hardware_set_soft_timer_rsp_t* gecko_cmd_hardware_set_soft_timer(uint32 time,uint8 handle,uint8 single_shot)
{
    
    gecko_cmd_msg->data.cmd_hardware_set_soft_timer.time=time;
    gecko_cmd_msg->data.cmd_hardware_set_soft_timer.handle=handle;
    gecko_cmd_msg->data.cmd_hardware_set_soft_timer.single_shot=single_shot;
    gecko_cmd_msg->header=(gecko_cmd_hardware_set_soft_timer_id+(((6)&0xff)<<8)+(((6)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_hardware_set_soft_timer;
}

/** 
*
* gecko_cmd_hardware_get_time
*
* Deprecated. Use Sleep Timer component (sl_sleeptimer.h) for the same
* functionality. Call sl_sleeptimer_get_tick_count64 to get current tick count.
* Sleep Timer provides APIs for conversions between ticks and milliseconds.
* 
* Get elapsed time since last reset. 
*
*
**/
BGLIB_DEPRECATED_API 
static inline struct gecko_msg_hardware_get_time_rsp_t* gecko_cmd_hardware_get_time()
{
    
    gecko_cmd_msg->header=(gecko_cmd_hardware_get_time_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_hardware_get_time;
}

/** 
*
* gecko_cmd_hardware_set_lazy_soft_timer
*
* Start a software timer with slack. The slack parameter allows the stack to
* optimize wakeups and save power. The timer event is triggered between time and
* time + slack. See also description of hardware_set_soft_timer command. 
*
* @param time   Interval between how often to send events in hardware clock ticks (1 second is
*  equal to 32768 ticks).
*  
*  The smallest interval value supported is 328, which is around 10 milliseconds.
*  Any parameters between 0 and 328 will be rounded up to 328. The maximum value
*  is 2147483647, which corresponds to about 18.2 hours.
*  
*  If time is 0, removes the scheduled timer with the same handle.
* @param slack   Slack time in hardware clock ticks
* @param handle   Timer handle to use, which is returned in timeout event
* @param single_shot   Timer mode. Values:
*  
*      0: false (timer is repeating)  
*       1: true (timer runs only once)
*
* Events generated
*
* gecko_evt_hardware_soft_timer - Sent after this timer has lapsed.
*
**/

static inline struct gecko_msg_hardware_set_lazy_soft_timer_rsp_t* gecko_cmd_hardware_set_lazy_soft_timer(uint32 time,uint32 slack,uint8 handle,uint8 single_shot)
{
    
    gecko_cmd_msg->data.cmd_hardware_set_lazy_soft_timer.time=time;
    gecko_cmd_msg->data.cmd_hardware_set_lazy_soft_timer.slack=slack;
    gecko_cmd_msg->data.cmd_hardware_set_lazy_soft_timer.handle=handle;
    gecko_cmd_msg->data.cmd_hardware_set_lazy_soft_timer.single_shot=single_shot;
    gecko_cmd_msg->header=(gecko_cmd_hardware_set_lazy_soft_timer_id+(((10)&0xff)<<8)+(((10)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_hardware_set_lazy_soft_timer;
}

/** 
*
* gecko_cmd_flash_ps_erase_all
*
* Delete all PS keys and their corresponding values. 
*
*
**/

static inline struct gecko_msg_flash_ps_erase_all_rsp_t* gecko_cmd_flash_ps_erase_all()
{
    
    gecko_cmd_msg->header=(gecko_cmd_flash_ps_erase_all_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_flash_ps_erase_all;
}

/** 
*
* gecko_cmd_flash_ps_save
*
* Store a value into the specified PS key. Allowed PS keys are in range from
* 0x4000 to 0x407F. At most, 56 bytes user data can be stored in one PS key.
* Error code 0x018a (command_too_long) is returned if the value data is more
* than 56 bytes. 
*
* @param key   PS key
* @param value_len   Array length
* @param value_data   Value to store into the specified PS key
*
**/

static inline struct gecko_msg_flash_ps_save_rsp_t* gecko_cmd_flash_ps_save(uint16 key,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 3)
    {
        gecko_rsp_msg->data.rsp_flash_ps_save.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_flash_ps_save;
    }

    
    gecko_cmd_msg->data.cmd_flash_ps_save.key=key;
    gecko_cmd_msg->data.cmd_flash_ps_save.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_flash_ps_save.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_flash_ps_save_id+(((3+value_len)&0xff)<<8)+(((3+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_flash_ps_save;
}

/** 
*
* gecko_cmd_flash_ps_load
*
* Retrieve the value of the specified PS key. 
*
* @param key   PS key of the value to be retrieved
*
**/

static inline struct gecko_msg_flash_ps_load_rsp_t* gecko_cmd_flash_ps_load(uint16 key)
{
    
    gecko_cmd_msg->data.cmd_flash_ps_load.key=key;
    gecko_cmd_msg->header=(gecko_cmd_flash_ps_load_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_flash_ps_load;
}

/** 
*
* gecko_cmd_flash_ps_erase
*
* Delete a single PS key and its value from the persistent store. 
*
* @param key   PS key to delete
*
**/

static inline struct gecko_msg_flash_ps_erase_rsp_t* gecko_cmd_flash_ps_erase(uint16 key)
{
    
    gecko_cmd_msg->data.cmd_flash_ps_erase.key=key;
    gecko_cmd_msg->header=(gecko_cmd_flash_ps_erase_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_flash_ps_erase;
}

/** 
*
* gecko_cmd_test_dtm_tx
*
* Start a transmitter test against a separate Bluetooth tester device. When the
* command is processed by the radio, a test_dtm_completed event is triggered.
* This event indicates whether the test started successfully.
* 
* In the transmitter test, the device sends packets continuously with a fixed
* interval. The type and length of each packet is set by packet_type and length
* parameters. The parameter phy specifies which PHY is used to transmit the
* packets. All devices support at least 1M PHY. A special packet type,
* test_pkt_carrier , can be used to transmit continuous unmodulated carrier. The
* length field is ignored in this mode.
* 
* The test may be stopped using the test_dtm_end command. 
*
* @param packet_type   Packet type to transmit
* @param length   Packet length in bytes
*  
*   Range: 0-255
* @param channel   Bluetooth channel
*  
*   Range: 0-39
*  
*  Channel is (F - 2402) / 2,
*  
*  where F is frequency in MHz
* @param phy   PHY to use
*
* Events generated
*
* gecko_evt_test_dtm_completed - This event is received when the command is processed.
*
**/

static inline struct gecko_msg_test_dtm_tx_rsp_t* gecko_cmd_test_dtm_tx(uint8 packet_type,uint8 length,uint8 channel,uint8 phy)
{
    
    gecko_cmd_msg->data.cmd_test_dtm_tx.packet_type=packet_type;
    gecko_cmd_msg->data.cmd_test_dtm_tx.length=length;
    gecko_cmd_msg->data.cmd_test_dtm_tx.channel=channel;
    gecko_cmd_msg->data.cmd_test_dtm_tx.phy=phy;
    gecko_cmd_msg->header=(gecko_cmd_test_dtm_tx_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_test_dtm_tx;
}

/** 
*
* gecko_cmd_test_dtm_rx
*
* Start a receiver test against a separate Bluetooth tester device. When the
* command is processed by the radio, a test_dtm_completed event is triggered.
* This event indicates whether the test started successfully.
* 
* Parameter phy specifies which PHY is used to receive the packets. All devices
* support at least 1M PHY.
* 
* The test may be stopped using the test_dtm_end command. This will trigger
* another test_dtm_completed event, which carries the number of packets received
* during the test. 
*
* @param channel   Bluetooth channel
*  
*   Range: 0-39
*  
*  Channel is (F - 2402) / 2,
*  
*  where F is frequency in MHz
* @param phy   PHY to use
*
* Events generated
*
* gecko_evt_test_dtm_completed - This event is received when the command is processed.
*
**/

static inline struct gecko_msg_test_dtm_rx_rsp_t* gecko_cmd_test_dtm_rx(uint8 channel,uint8 phy)
{
    
    gecko_cmd_msg->data.cmd_test_dtm_rx.channel=channel;
    gecko_cmd_msg->data.cmd_test_dtm_rx.phy=phy;
    gecko_cmd_msg->header=(gecko_cmd_test_dtm_rx_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_test_dtm_rx;
}

/** 
*
* gecko_cmd_test_dtm_end
*
* End a transmitter or a receiver test. When the command is processed by the
* radio and the test has ended, a test_dtm_completed event is triggered. 
*
*
* Events generated
*
* gecko_evt_test_dtm_completed - Received when the command is processed by the radio and the test has ended.
*
**/

static inline struct gecko_msg_test_dtm_end_rsp_t* gecko_cmd_test_dtm_end()
{
    
    gecko_cmd_msg->header=(gecko_cmd_test_dtm_end_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_test_dtm_end;
}

/** 
*
* gecko_cmd_sm_set_bondable_mode
*
* Set whether the device should accept new bondings. By default, the device does
* not accept new bondings. 
*
* @param bondable   Bondable mode. Values:
*  
*      0: New bondings not accepted  
*      1: Bondings allowed  
*  
*  Default value: 0
*
**/

static inline struct gecko_msg_sm_set_bondable_mode_rsp_t* gecko_cmd_sm_set_bondable_mode(uint8 bondable)
{
    
    gecko_cmd_msg->data.cmd_sm_set_bondable_mode.bondable=bondable;
    gecko_cmd_msg->header=(gecko_cmd_sm_set_bondable_mode_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_set_bondable_mode;
}

/** 
*
* gecko_cmd_sm_configure
*
* Configure security requirements and I/O capabilities of the system. 
*
* @param flags   Security requirement bitmask.
*  
*  Bit 0:
*  
*      0: Allow bonding without MITM protection  
*      1: Bonding requires MITM protection  
*  
*  Bit 1:
*  
*      0: Allow encryption without bonding  
*      1: Encryption requires bonding. Note that this setting will also enable
*      bonding.  
*  
*  Bit 2:
*  
*      0: Allow bonding with legacy pairing  
*      1: Secure connections only  
*  
*  Bit 3:
*  
*      0: Bonding request does not need to be confirmed  
*      1: Bonding requests need to be confirmed. Received bonding requests are
*      notified by sm_confirm_bonding events.  
*  
*  Bit 4:
*  
*      0: Allow all connections  
*      1: Allow connections only from bonded devices  
*  
*  Bit 5 to 7: Reserved
*  
*  Default value: 0x00
* @param io_capabilities   I/O Capabilities. See link.
*
**/

static inline struct gecko_msg_sm_configure_rsp_t* gecko_cmd_sm_configure(uint8 flags,uint8 io_capabilities)
{
    
    gecko_cmd_msg->data.cmd_sm_configure.flags=flags;
    gecko_cmd_msg->data.cmd_sm_configure.io_capabilities=io_capabilities;
    gecko_cmd_msg->header=(gecko_cmd_sm_configure_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_configure;
}

/** 
*
* gecko_cmd_sm_store_bonding_configuration
*
* Set the maximum allowed bonding count and bonding policy. The maximum number
* of bondings that can be supported depends on how much user data is stored in
* the NVM and the NVM size. When bond policy value 1 or 2 is selected the stack
* will automatically write the new bond, as per the policy, only if the maximum
* allowed bonding count has been reached. If the stack is not able to write a
* new bond for any other reason (e.g. nvm full) then an error will be thrown
* through the bonding_failed event indicating why the bonding could not be
* written. It is left up to the application to manually release space from the
* nvm (e.g. by deleting one of the existing bonds or application data) so that a
* new bond can be saved. The default value is 13. 
*
* @param max_bonding_count   Maximum allowed bonding count. Range: 1 to 32
* @param policy_flags   Bonding policy. Values:
*  
*      0: If database is full, new bonding attempts will fail  
*      1: New bonding will overwrite the oldest existing bonding  
*      2: New bonding will overwrite the existing bonding that was used the
*      longest time ago  
*  
*  Default: 0
*
**/

static inline struct gecko_msg_sm_store_bonding_configuration_rsp_t* gecko_cmd_sm_store_bonding_configuration(uint8 max_bonding_count,uint8 policy_flags)
{
    
    gecko_cmd_msg->data.cmd_sm_store_bonding_configuration.max_bonding_count=max_bonding_count;
    gecko_cmd_msg->data.cmd_sm_store_bonding_configuration.policy_flags=policy_flags;
    gecko_cmd_msg->header=(gecko_cmd_sm_store_bonding_configuration_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_store_bonding_configuration;
}

/** 
*
* gecko_cmd_sm_increase_security
*
* Enhance the security of a connection to current security requirements. On an
* unencrypted connection, it will encrypt the connection and will also perform
* bonding if requested by both devices. On an encrypted connection, it will
* cause the connection to be re-encrypted. 
*
* @param connection   Connection handle
*
* Events generated
*
* gecko_evt_le_connection_parameters - Triggered after increasing security has been completed successfully and
*  indicates the latest security mode of the connection.
* gecko_evt_sm_bonded - Triggered if pairing or bonding was performed in this operation and the result
*  is successful.
* gecko_evt_sm_bonding_failed - Triggered if pairing or bonding was performed in this operation and the result
*  has failed.
*
**/

static inline struct gecko_msg_sm_increase_security_rsp_t* gecko_cmd_sm_increase_security(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_sm_increase_security.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_sm_increase_security_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_increase_security;
}

/** 
*
* gecko_cmd_sm_delete_bonding
*
* Delete specified bonding information or whitelist from the persistent store. 
*
* @param bonding   Bonding handle
*
**/

static inline struct gecko_msg_sm_delete_bonding_rsp_t* gecko_cmd_sm_delete_bonding(uint8 bonding)
{
    
    gecko_cmd_msg->data.cmd_sm_delete_bonding.bonding=bonding;
    gecko_cmd_msg->header=(gecko_cmd_sm_delete_bonding_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_delete_bonding;
}

/** 
*
* gecko_cmd_sm_delete_bondings
*
* Delete all bonding information and whitelist from the persistent store. 
*
*
**/

static inline struct gecko_msg_sm_delete_bondings_rsp_t* gecko_cmd_sm_delete_bondings()
{
    
    gecko_cmd_msg->header=(gecko_cmd_sm_delete_bondings_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_delete_bondings;
}

/** 
*
* gecko_cmd_sm_enter_passkey
*
* Enter a passkey after receiving a passkey request event. 
*
* @param connection   Connection handle
* @param passkey   Passkey. Valid range: 0-999999. Set -1 to cancel pairing.
*
**/

static inline struct gecko_msg_sm_enter_passkey_rsp_t* gecko_cmd_sm_enter_passkey(uint8 connection,int32 passkey)
{
    
    gecko_cmd_msg->data.cmd_sm_enter_passkey.connection=connection;
    gecko_cmd_msg->data.cmd_sm_enter_passkey.passkey=passkey;
    gecko_cmd_msg->header=(gecko_cmd_sm_enter_passkey_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_enter_passkey;
}

/** 
*
* gecko_cmd_sm_passkey_confirm
*
* Accept or reject the reported passkey confirm value. 
*
* @param connection   Connection handle
* @param confirm   Acceptance. Values:
*  
*       0: Reject  
*       1: Accept confirm value
*
**/

static inline struct gecko_msg_sm_passkey_confirm_rsp_t* gecko_cmd_sm_passkey_confirm(uint8 connection,uint8 confirm)
{
    
    gecko_cmd_msg->data.cmd_sm_passkey_confirm.connection=connection;
    gecko_cmd_msg->data.cmd_sm_passkey_confirm.confirm=confirm;
    gecko_cmd_msg->header=(gecko_cmd_sm_passkey_confirm_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_passkey_confirm;
}

/** 
*
* gecko_cmd_sm_set_oob_data
*
* Set OOB data (out-of-band encryption data) for legacy pairing for a device.
* OOB data may be, for example, a PIN code exchanged over an alternate path,
* such as NFC. The device will not allow any other bonding if OOB data is set.
* OOB data can't be set simultaneously with secure connections OOB data. 
*
* @param oob_data_len   Array length
* @param oob_data_data   OOB data. To set OOB data, send a 16-byte array. To clear OOB data, send a
*  zero-length array.
*
**/

static inline struct gecko_msg_sm_set_oob_data_rsp_t* gecko_cmd_sm_set_oob_data(uint8 oob_data_len, const uint8* oob_data_data)
{
    if ((uint16_t)oob_data_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_sm_set_oob_data.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_sm_set_oob_data;
    }

    
    gecko_cmd_msg->data.cmd_sm_set_oob_data.oob_data.len=oob_data_len;
    memcpy(gecko_cmd_msg->data.cmd_sm_set_oob_data.oob_data.data,oob_data_data,oob_data_len);
    gecko_cmd_msg->header=(gecko_cmd_sm_set_oob_data_id+(((1+oob_data_len)&0xff)<<8)+(((1+oob_data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_set_oob_data;
}

/** 
*
* gecko_cmd_sm_list_all_bondings
*
* List all bondings stored in the bonding database. Bondings are reported by the
* sm_list_bonding_entry event for each bonding and the report is ended with
* sm_list_all_bondings_complete event. Use only for debugging purposes because
* reading from the persistent store is relatively slow. 
*
*
* Events generated
*
* gecko_evt_sm_list_bonding_entry - 
* gecko_evt_sm_list_all_bondings_complete - 
*
**/

static inline struct gecko_msg_sm_list_all_bondings_rsp_t* gecko_cmd_sm_list_all_bondings()
{
    
    gecko_cmd_msg->header=(gecko_cmd_sm_list_all_bondings_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_list_all_bondings;
}

/** 
*
* gecko_cmd_sm_bonding_confirm
*
* Accept or reject the bonding request. 
*
* @param connection   Connection handle
* @param confirm   Acceptance. Values:
*  
*       0: Reject  
*       1: Accept bonding request
*
**/

static inline struct gecko_msg_sm_bonding_confirm_rsp_t* gecko_cmd_sm_bonding_confirm(uint8 connection,uint8 confirm)
{
    
    gecko_cmd_msg->data.cmd_sm_bonding_confirm.connection=connection;
    gecko_cmd_msg->data.cmd_sm_bonding_confirm.confirm=confirm;
    gecko_cmd_msg->header=(gecko_cmd_sm_bonding_confirm_id+(((2)&0xff)<<8)+(((2)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_bonding_confirm;
}

/** 
*
* gecko_cmd_sm_set_debug_mode
*
* Set Security Manager in debug mode. In this mode, the secure connections
* bonding uses known debug keys, so that the encrypted packet can be opened by
* Bluetooth protocol analyzer. To disable the debug mode, restart the device.
* 
* Bondings made in debug mode are unsecure. 
*
*
**/

static inline struct gecko_msg_sm_set_debug_mode_rsp_t* gecko_cmd_sm_set_debug_mode()
{
    
    gecko_cmd_msg->header=(gecko_cmd_sm_set_debug_mode_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_set_debug_mode;
}

/** 
*
* gecko_cmd_sm_set_passkey
*
* Enter a fixed passkey, which will be used in the sm_passkey_display event. 
*
* @param passkey   Passkey. Valid range: 0-999999. Set -1 to disable and start using random
*  passkeys.
*
**/

static inline struct gecko_msg_sm_set_passkey_rsp_t* gecko_cmd_sm_set_passkey(int32 passkey)
{
    
    gecko_cmd_msg->data.cmd_sm_set_passkey.passkey=passkey;
    gecko_cmd_msg->header=(gecko_cmd_sm_set_passkey_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_set_passkey;
}

/** 
*
* gecko_cmd_sm_use_sc_oob
*
* Enable the use of OOB data (out-of-band encryption data) for a device for
* secure connections pairing. Enabling will generate new OOB data and confirm
* values, which can be sent to the remote device. After enabling the secure
* connections OOB data, the remote devices OOB data can be set with
* sm_set_sc_remote_oob_data. Calling this function will erase any set remote
* device OOB data and confirm values. The device will not allow any other
* bonding if OOB data is set. The secure connections OOB data cannot be enabled
* simultaneously with legacy pairing OOB data. 
*
* @param enable   Enable OOB with secure connections pairing. Values:
*  
*      0: disable  
*       1: enable
*
**/

static inline struct gecko_msg_sm_use_sc_oob_rsp_t* gecko_cmd_sm_use_sc_oob(uint8 enable)
{
    
    gecko_cmd_msg->data.cmd_sm_use_sc_oob.enable=enable;
    gecko_cmd_msg->header=(gecko_cmd_sm_use_sc_oob_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_use_sc_oob;
}

/** 
*
* gecko_cmd_sm_set_sc_remote_oob_data
*
* Set OOB data and confirm values (out-of-band encryption) received from the
* remote device for secure connections pairing. OOB data must be enabled with
* sm_use_sc_oob before setting the remote device OOB data. 
*
* @param oob_data_len   Array length
* @param oob_data_data   Remote device OOB data and confirm values. To set OOB data, send a 32-byte
*  array. First 16-bytes is OOB data and last 16-bytes the confirm value. To
*  clear OOB data, send a zero-length array.
*
**/

static inline struct gecko_msg_sm_set_sc_remote_oob_data_rsp_t* gecko_cmd_sm_set_sc_remote_oob_data(uint8 oob_data_len, const uint8* oob_data_data)
{
    if ((uint16_t)oob_data_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_sm_set_sc_remote_oob_data.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_sm_set_sc_remote_oob_data;
    }

    
    gecko_cmd_msg->data.cmd_sm_set_sc_remote_oob_data.oob_data.len=oob_data_len;
    memcpy(gecko_cmd_msg->data.cmd_sm_set_sc_remote_oob_data.oob_data.data,oob_data_data,oob_data_len);
    gecko_cmd_msg->header=(gecko_cmd_sm_set_sc_remote_oob_data_id+(((1+oob_data_len)&0xff)<<8)+(((1+oob_data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_set_sc_remote_oob_data;
}

/** 
*
* gecko_cmd_sm_add_to_whitelist
*
* Add device to whitelist, which can be enabled with le_gap_enable_whitelisting. 
*
* @param address   Address of the device added to whitelist
* @param address_type   Address type of the device added to whitelist
*
**/

static inline struct gecko_msg_sm_add_to_whitelist_rsp_t* gecko_cmd_sm_add_to_whitelist(bd_addr address,uint8 address_type)
{
    
    memcpy(&gecko_cmd_msg->data.cmd_sm_add_to_whitelist.address,&address,sizeof(bd_addr));
    gecko_cmd_msg->data.cmd_sm_add_to_whitelist.address_type=address_type;
    gecko_cmd_msg->header=(gecko_cmd_sm_add_to_whitelist_id+(((7)&0xff)<<8)+(((7)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_add_to_whitelist;
}

/** 
*
* gecko_cmd_sm_set_minimum_key_size
*
* Set the minimum allowed key size used for bonding. The default value is 16
* bytes. 
*
* @param minimum_key_size   Minimum allowed key size for bonding. Range: 7 to 16
*
**/

static inline struct gecko_msg_sm_set_minimum_key_size_rsp_t* gecko_cmd_sm_set_minimum_key_size(uint8 minimum_key_size)
{
    
    gecko_cmd_msg->data.cmd_sm_set_minimum_key_size.minimum_key_size=minimum_key_size;
    gecko_cmd_msg->header=(gecko_cmd_sm_set_minimum_key_size_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_sm_set_minimum_key_size;
}

/** 
*
* gecko_cmd_homekit_configure
*
* Configure the Apple HomeKit accessory and its settings. The configuration can
* be reinitialized at run time. New fast advertising parameters will be used for
* next fast advertising. 
*
* @param i2c_address   I2C address of Apple authentication coprocessor
* @param support_display   A flag to indicate that the display support is enabled in the accessory.  
*  A pin code will be randomly generated during the pairing process and
*  homekit_setupcode_display event will be triggered to ensure that the pin code
*  can be displayed.
*  
*      0: Display support disabled  
*      1: Display support enabled
* @param hap_attribute_features   The value of Apple HomeKit pairing features supported in pairing service
*  feature characteristic.
*  
*      0x01: Supports Apple Authentication Coprocessor  
*      0x02: Supports Software Authentication  
*      0x00: Only for testing purposes when authentication methods are not
*      available. The accessory will be discovered as non-authenticated  
*      other: Reserved
* @param category   Apple HomeKit accessory category
* @param configuration_number   Apple HomeKit configuration number  
*  By default, this starts at 1. Accessories must increment this value after a
*  firmware update. This value must be managed by the application.
* @param fast_advert_interval   Fast advertising interval.  
*  The interval is used during fast advertising in disconnected state after
*  calling command homekit_event_notification when broadcast events advertising
*  is complete.
* @param fast_advert_timeout   Fast advertising timeout.  
*  The timeout is used during fast advertising in disconnected state after
*  calling command homekit_event_notification when broadcast events advertising
*  is complete.
*  
*      Time = Value x 100 ms
* @param flag   Apple HomeKit library configuration flag.
*  
*      0x00000001: Manual Bluetooth disconnection in HomeKit error case. When
*      enabled, a homekit_disconnection_required event will be produced when a
*      HomeKit error occurs.  
*      0x00000002: Manual set of scan response data. When enabled, use
*      le_gap_bt5_set_adv_data command to set custom scan response data. Also,
*      HomeKit library uses it to set the accessory local name.  
*      other: Reserved. Must be 0.
* @param broadcast_advert_timeout   Broadcast events advertising timeout.  
*  The timeout is used during broadcast events advertising in disconnected state
*  after calling command homekit_event_notification
*  
*      Time = Value x 100 ms
* @param model_name_len   Array length
* @param model_name_data   Model Name characteristic value from HomeKit Accessory Information service.
*  Mandatory for HomeKit software authentication usage.
*
**/

static inline struct gecko_msg_homekit_configure_rsp_t* gecko_cmd_homekit_configure(uint8 i2c_address,uint8 support_display,uint8 hap_attribute_features,uint16 category,uint8 configuration_number,uint16 fast_advert_interval,uint16 fast_advert_timeout,uint32 flag,uint16 broadcast_advert_timeout,uint8 model_name_len, const uint8* model_name_data)
{
    if ((uint16_t)model_name_len > BGLIB_MSG_MAX_PAYLOAD - 17)
    {
        gecko_rsp_msg->data.rsp_homekit_configure.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_homekit_configure;
    }

    
    gecko_cmd_msg->data.cmd_homekit_configure.i2c_address=i2c_address;
    gecko_cmd_msg->data.cmd_homekit_configure.support_display=support_display;
    gecko_cmd_msg->data.cmd_homekit_configure.hap_attribute_features=hap_attribute_features;
    gecko_cmd_msg->data.cmd_homekit_configure.category=category;
    gecko_cmd_msg->data.cmd_homekit_configure.configuration_number=configuration_number;
    gecko_cmd_msg->data.cmd_homekit_configure.fast_advert_interval=fast_advert_interval;
    gecko_cmd_msg->data.cmd_homekit_configure.fast_advert_timeout=fast_advert_timeout;
    gecko_cmd_msg->data.cmd_homekit_configure.flag=flag;
    gecko_cmd_msg->data.cmd_homekit_configure.broadcast_advert_timeout=broadcast_advert_timeout;
    gecko_cmd_msg->data.cmd_homekit_configure.model_name.len=model_name_len;
    memcpy(gecko_cmd_msg->data.cmd_homekit_configure.model_name.data,model_name_data,model_name_len);
    gecko_cmd_msg->header=(gecko_cmd_homekit_configure_id+(((17+model_name_len)&0xff)<<8)+(((17+model_name_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_configure;
}

/** 
*
* gecko_cmd_homekit_advertise
*
* Start or stop Apple HomeKit accessory advertising. The command and parameters
* will take effect immediately. If the given parameters can't be used in the
* currently active mode, an error is returned. 
*
* @param enable   A flag to enable or disable Apple HomeKit advertising
*  
*      1: Enable advertising  
*      0: Disable advertising
* @param interval_min   Minimum advertising interval. See GAP command le_gap_set_adv_parameters
* @param interval_max   Maximum advertising interval. See GAP command le_gap_set_adv_parameters
* @param channel_map   Advertising channel map. See GAP command le_gap_set_adv_parameters
*
**/

static inline struct gecko_msg_homekit_advertise_rsp_t* gecko_cmd_homekit_advertise(uint8 enable,uint16 interval_min,uint16 interval_max,uint8 channel_map)
{
    
    gecko_cmd_msg->data.cmd_homekit_advertise.enable=enable;
    gecko_cmd_msg->data.cmd_homekit_advertise.interval_min=interval_min;
    gecko_cmd_msg->data.cmd_homekit_advertise.interval_max=interval_max;
    gecko_cmd_msg->data.cmd_homekit_advertise.channel_map=channel_map;
    gecko_cmd_msg->header=(gecko_cmd_homekit_advertise_id+(((6)&0xff)<<8)+(((6)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_advertise;
}

/** 
*
* gecko_cmd_homekit_delete_pairings
*
* Delete all Apple HomeKit pairing data. Additionally, it resets all required
* HomeKit settings to factory state, e.g., it resets GSN value and generates a
* new Device ID. 
*
*
**/

static inline struct gecko_msg_homekit_delete_pairings_rsp_t* gecko_cmd_homekit_delete_pairings()
{
    
    gecko_cmd_msg->header=(gecko_cmd_homekit_delete_pairings_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_delete_pairings;
}

/** 
*
* gecko_cmd_homekit_check_authcp
*
* Create an I2C test connection with Apple authentication co-processor and
* return an error if communication fails. 
*
*
**/

static inline struct gecko_msg_homekit_check_authcp_rsp_t* gecko_cmd_homekit_check_authcp()
{
    
    gecko_cmd_msg->header=(gecko_cmd_homekit_check_authcp_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_check_authcp;
}

/** 
*
* gecko_cmd_homekit_get_pairing_id
*
* Get pairing ID of the connected iOS device. 
*
* @param connection   
*
**/

static inline struct gecko_msg_homekit_get_pairing_id_rsp_t* gecko_cmd_homekit_get_pairing_id(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_homekit_get_pairing_id.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_homekit_get_pairing_id_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_get_pairing_id;
}

/** 
*
* gecko_cmd_homekit_send_write_response
*
* Send a response to a homekit_write_request event. The response needs to be
* sent within 30 seconds, otherwise additional GATT transactions are not allowed
* by the remote side.  
*   
* If the status_code is set to 0, the HAP will send a response informing that
* the write operation was processed successfully and other values will cause the
* HAP to send a HAP error status response. 
*
* @param connection   
* @param characteristic   
* @param status_code   HomeKit status code.
*
**/

static inline struct gecko_msg_homekit_send_write_response_rsp_t* gecko_cmd_homekit_send_write_response(uint8 connection,uint16 characteristic,uint8 status_code)
{
    
    gecko_cmd_msg->data.cmd_homekit_send_write_response.connection=connection;
    gecko_cmd_msg->data.cmd_homekit_send_write_response.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_homekit_send_write_response.status_code=status_code;
    gecko_cmd_msg->header=(gecko_cmd_homekit_send_write_response_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_send_write_response;
}

/** 
*
* gecko_cmd_homekit_send_read_response
*
* Send a response to a homekit_read_request event. The response needs to be sent
* within 30 seconds, otherwise further GATT transactions are not allowed by the
* remote side.  
*   
* If status_code is set to 0, the characteristic value is sent to the remote
* GATT client through HomeKit library in a standard way. Other status_code
* values cause a HAP error status response instead of sending data.  
*   
* If the value data size is less than attribute_size, the Apple HomeKit library
* will send new homekit_read_request event with a suitable offset. The Apple
* HomeKit library provides automatic formatting for both the frame and
* encryption. 
*
* @param connection   
* @param characteristic   
* @param status_code   HomeKit Status Code
* @param attribute_size   Size of attribute value
* @param value_len   Array length
* @param value_data   Characteristic value to send to the GATT client through the HomeKit library.
*  This is ignored if status_code is not set to 0.
*
**/

static inline struct gecko_msg_homekit_send_read_response_rsp_t* gecko_cmd_homekit_send_read_response(uint8 connection,uint16 characteristic,uint8 status_code,uint16 attribute_size,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 7)
    {
        gecko_rsp_msg->data.rsp_homekit_send_read_response.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_homekit_send_read_response;
    }

    
    gecko_cmd_msg->data.cmd_homekit_send_read_response.connection=connection;
    gecko_cmd_msg->data.cmd_homekit_send_read_response.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_homekit_send_read_response.status_code=status_code;
    gecko_cmd_msg->data.cmd_homekit_send_read_response.attribute_size=attribute_size;
    gecko_cmd_msg->data.cmd_homekit_send_read_response.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_homekit_send_read_response.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_homekit_send_read_response_id+(((7+value_len)&0xff)<<8)+(((7+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_send_read_response;
}

/** 
*
* gecko_cmd_homekit_gsn_action
*
* Reset or store the GSN (Global State Number) value. 
*
* @param action   Actions:
*  
*      0: Reset GSN value to default state  
*      1: Store GSN value to a PS-key (flash)
*
**/

static inline struct gecko_msg_homekit_gsn_action_rsp_t* gecko_cmd_homekit_gsn_action(uint8 action)
{
    
    gecko_cmd_msg->data.cmd_homekit_gsn_action.action=action;
    gecko_cmd_msg->header=(gecko_cmd_homekit_gsn_action_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_gsn_action;
}

/** 
*
* gecko_cmd_homekit_event_notification
*
* Take an appropriate action according to connection state and the originator of
* the change. When device is connected and a local change occurs, it sends an
* empty indication to the controller. When device is disconnected, it starts
* broadcast events advertising. After timeout, it starts fast advertising.
* Broadcast and fast advertising parameters are set in homekit_configure. After
* fast advertising timeout, it reverts to previous advertising settings. For
* both states, it sets the appropriate Global State Number value according to
* HomeKit specification rules. 
*
* @param connection   Connection handle. Ignored for disconnected state.
* @param characteristic   
* @param change_originator   Origin of the characteristic value change:
*  
*      0: Remote change (from controller)  
*      1: Local change (from accessory)
* @param value_len   Array length
* @param value_data   Broadcast notify value.
*
**/

static inline struct gecko_msg_homekit_event_notification_rsp_t* gecko_cmd_homekit_event_notification(uint8 connection,uint16 characteristic,uint8 change_originator,uint8 value_len, const uint8* value_data)
{
    if ((uint16_t)value_len > BGLIB_MSG_MAX_PAYLOAD - 5)
    {
        gecko_rsp_msg->data.rsp_homekit_event_notification.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_homekit_event_notification;
    }

    
    gecko_cmd_msg->data.cmd_homekit_event_notification.connection=connection;
    gecko_cmd_msg->data.cmd_homekit_event_notification.characteristic=characteristic;
    gecko_cmd_msg->data.cmd_homekit_event_notification.change_originator=change_originator;
    gecko_cmd_msg->data.cmd_homekit_event_notification.value.len=value_len;
    memcpy(gecko_cmd_msg->data.cmd_homekit_event_notification.value.data,value_data,value_len);
    gecko_cmd_msg->header=(gecko_cmd_homekit_event_notification_id+(((5+value_len)&0xff)<<8)+(((5+value_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_event_notification;
}

/** 
*
* gecko_cmd_homekit_broadcast_action
*
* Delete or store broadcast advertising data, as shown below. 
*
* @param action   Actions:
*  
*      0x00: Delete broadcast advertising data. No additional parameters are
*      required  
*      0x01: Store broadcast advertising data (key, characteristics
*      configuration) to non volatile memory. No additional parameters are
*      required  
*      other: Reserved
* @param params_len   Array length
* @param params_data   Additional parameters for action.
*
**/

static inline struct gecko_msg_homekit_broadcast_action_rsp_t* gecko_cmd_homekit_broadcast_action(uint8 action,uint8 params_len, const uint8* params_data)
{
    if ((uint16_t)params_len > BGLIB_MSG_MAX_PAYLOAD - 2)
    {
        gecko_rsp_msg->data.rsp_homekit_broadcast_action.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_homekit_broadcast_action;
    }

    
    gecko_cmd_msg->data.cmd_homekit_broadcast_action.action=action;
    gecko_cmd_msg->data.cmd_homekit_broadcast_action.params.len=params_len;
    memcpy(gecko_cmd_msg->data.cmd_homekit_broadcast_action.params.data,params_data,params_len);
    gecko_cmd_msg->header=(gecko_cmd_homekit_broadcast_action_id+(((2+params_len)&0xff)<<8)+(((2+params_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_broadcast_action;
}

/** 
*
* gecko_cmd_homekit_configure_product_data
*
* Configure the Apple HomeKit library. This is additional configuration
* introduced in HAP revision R15. 
*
* @param product_data_len   Array length
* @param product_data_data   Product Data characteristic value from HomeKit Accessory Information service.
*  Mandatory (HAP revision R15 and later).
*
**/

static inline struct gecko_msg_homekit_configure_product_data_rsp_t* gecko_cmd_homekit_configure_product_data(uint8 product_data_len, const uint8* product_data_data)
{
    if ((uint16_t)product_data_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_homekit_configure_product_data.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_homekit_configure_product_data;
    }

    
    gecko_cmd_msg->data.cmd_homekit_configure_product_data.product_data.len=product_data_len;
    memcpy(gecko_cmd_msg->data.cmd_homekit_configure_product_data.product_data.data,product_data_data,product_data_len);
    gecko_cmd_msg->header=(gecko_cmd_homekit_configure_product_data_id+(((1+product_data_len)&0xff)<<8)+(((1+product_data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_homekit_configure_product_data;
}

/** 
*
* gecko_cmd_coex_set_options
*
* Configure coexistence options at runtime. 
*
* @param mask   Mask defines which coexistence options are changed.
* @param options   Value of options to be changed. This parameter is used together with the mask
*  parameter.
*
**/

static inline struct gecko_msg_coex_set_options_rsp_t* gecko_cmd_coex_set_options(uint32 mask,uint32 options)
{
    
    gecko_cmd_msg->data.cmd_coex_set_options.mask=mask;
    gecko_cmd_msg->data.cmd_coex_set_options.options=options;
    gecko_cmd_msg->header=(gecko_cmd_coex_set_options_id+(((8)&0xff)<<8)+(((8)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_coex_set_options;
}

/** 
*
* gecko_cmd_coex_get_counters
*
* Read coexistence statistic counters from the device. Response contains the
* list of uint32 type counter values. Counters in the list are in following
* order: low priority requested, high priority requested, low priority denied,
* high priority denied, low-priority TX aborted, and high-priority TX aborted.
* Passing a non-zero value also resets counters. 
*
* @param reset   Reset counters if parameter value is not zero.
*
**/

static inline struct gecko_msg_coex_get_counters_rsp_t* gecko_cmd_coex_get_counters(uint8 reset)
{
    
    gecko_cmd_msg->data.cmd_coex_get_counters.reset=reset;
    gecko_cmd_msg->header=(gecko_cmd_coex_get_counters_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_coex_get_counters;
}

/** 
*
* gecko_cmd_coex_set_parameters
*
* Configure coexistence parameters. 
*
* @param priority   Coexistence priority threshold. Coexistence priority is toggled if priority is
*  below this value.
* @param request   Coexistence request threshold. Coexistence request is toggled if priority is
*  below this value.
* @param pwm_period   PWM functionality period length in 1ms units
* @param pwm_dutycycle   PWM functionality dutycycle in percentage
*
**/

static inline struct gecko_msg_coex_set_parameters_rsp_t* gecko_cmd_coex_set_parameters(uint8 priority,uint8 request,uint8 pwm_period,uint8 pwm_dutycycle)
{
    
    gecko_cmd_msg->data.cmd_coex_set_parameters.priority=priority;
    gecko_cmd_msg->data.cmd_coex_set_parameters.request=request;
    gecko_cmd_msg->data.cmd_coex_set_parameters.pwm_period=pwm_period;
    gecko_cmd_msg->data.cmd_coex_set_parameters.pwm_dutycycle=pwm_dutycycle;
    gecko_cmd_msg->header=(gecko_cmd_coex_set_parameters_id+(((4)&0xff)<<8)+(((4)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_coex_set_parameters;
}

/** 
*
* gecko_cmd_coex_set_directional_priority_pulse
*
* Set Directional Priority Pulse Width 
*
* @param pulse   Directional priority pulse width in us
*
**/

static inline struct gecko_msg_coex_set_directional_priority_pulse_rsp_t* gecko_cmd_coex_set_directional_priority_pulse(uint8 pulse)
{
    
    gecko_cmd_msg->data.cmd_coex_set_directional_priority_pulse.pulse=pulse;
    gecko_cmd_msg->header=(gecko_cmd_coex_set_directional_priority_pulse_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_coex_set_directional_priority_pulse;
}

/** 
*
* gecko_cmd_l2cap_coc_send_connection_request
*
* Send LE credit-based connection request. 
*
* @param connection   Handle of the LE connection to be used for opening connection-oriented
*  channel.
* @param le_psm   LE Protocol/Service Multiplexer - LE_PSM
* @param mtu   The maximum size of payload data that the application on the device sending
*  the request can accept, i.e., the MTU corresponds to the maximum SDU size.
*  
*  Range: 23 to 65533.
*  
*  Application needs to handle segmentation and reassembly from PDU to SDU.
* @param mps   The maximum size of payload data that the L2CAP layer on the device sending
*  the request can accept, i.e., the MPS corresponds to the maximum PDU payload
*  size.
*  
*  Range: 23 to 250.
*  
*  That is the maximum size of data that the application can send using
*  l2cap_coc_send_data command or receive by l2cap_coc_data event.
* @param initial_credit   The initial credit value indicates the number of PDUs that the peer device can
*  send.
*
* Events generated
*
* gecko_evt_l2cap_coc_connection_response - Triggered when a LE credit-based connection connection response has been
*  received in response to this command.
* gecko_evt_l2cap_coc_channel_disconnected - Triggered when a LE credit-based connection connection response has not been
*  received within the 30 seconds timeout in response to this command.
*
**/

static inline struct gecko_msg_l2cap_coc_send_connection_request_rsp_t* gecko_cmd_l2cap_coc_send_connection_request(uint8 connection,uint16 le_psm,uint16 mtu,uint16 mps,uint16 initial_credit)
{
    
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_request.connection=connection;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_request.le_psm=le_psm;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_request.mtu=mtu;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_request.mps=mps;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_request.initial_credit=initial_credit;
    gecko_cmd_msg->header=(gecko_cmd_l2cap_coc_send_connection_request_id+(((9)&0xff)<<8)+(((9)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_l2cap_coc_send_connection_request;
}

/** 
*
* gecko_cmd_l2cap_coc_send_connection_response
*
* Send LE credit-based connection response. 
*
* @param connection   Handle of the LE connection to be used for opening connection-oriented channel
* @param cid   The CID represents the destination channel endpoint of the device sending the
*  response which is same as source CID field of corresponding request message
* @param mtu   The maximum size of payload data that the application on the device sending
*  the response can accept, i.e., the MTU corresponds to the maximum SDU size.
*  
*  Range: 23 to 65533.
*  
*  Application needs to handle segmentation and reassembly from PDU to SDU.
* @param mps   The maximum size of payload data that the L2CAP layer on the device sending
*  the response can accept, i.e., the MPS corresponds to the maximum PDU payload
*  size.
*  
*  Range: 23 to 250.
*  
*  That is the maximum size of data that the application is able to send using
*  l2cap_coc_send_data command or receive by l2cap_coc_data event.
* @param initial_credit   The initial credit value indicates the number of PDUs that the peer device can
*  send
* @param result   The result field indicates the outcome of the connection request.
*
**/

static inline struct gecko_msg_l2cap_coc_send_connection_response_rsp_t* gecko_cmd_l2cap_coc_send_connection_response(uint8 connection,uint16 cid,uint16 mtu,uint16 mps,uint16 initial_credit,uint16 result)
{
    
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_response.connection=connection;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_response.cid=cid;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_response.mtu=mtu;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_response.mps=mps;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_response.initial_credit=initial_credit;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_connection_response.result=result;
    gecko_cmd_msg->header=(gecko_cmd_l2cap_coc_send_connection_response_id+(((11)&0xff)<<8)+(((11)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_l2cap_coc_send_connection_response;
}

/** 
*
* gecko_cmd_l2cap_coc_send_le_flow_control_credit
*
* Send LE flow control credit indicating that the channel endpoint on local
* device is capable of receiving more data. 
*
* @param connection   Handle of the LE connection for sending flow control credit.
* @param cid   The CID represents the destination channel endpoint of the device sending the
*  flow control credit.
* @param credits   The credit value indicates the additional number of PDUs that the peer device
*  can send.
*
**/

static inline struct gecko_msg_l2cap_coc_send_le_flow_control_credit_rsp_t* gecko_cmd_l2cap_coc_send_le_flow_control_credit(uint8 connection,uint16 cid,uint16 credits)
{
    
    gecko_cmd_msg->data.cmd_l2cap_coc_send_le_flow_control_credit.connection=connection;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_le_flow_control_credit.cid=cid;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_le_flow_control_credit.credits=credits;
    gecko_cmd_msg->header=(gecko_cmd_l2cap_coc_send_le_flow_control_credit_id+(((5)&0xff)<<8)+(((5)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_l2cap_coc_send_le_flow_control_credit;
}

/** 
*
* gecko_cmd_l2cap_coc_send_disconnection_request
*
* Send L2CAP connection-oriented channel disconnection request. 
*
* @param connection   Handle of the LE connection for terminating the connection-oriented channel
* @param cid   The CID represents the destination channel endpoint of the device sending the
*  disconnection request.
*
* Events generated
*
* gecko_evt_l2cap_coc_channel_disconnected - Triggered when a L2CAP channel is disconnected in response to this command.
*
**/

static inline struct gecko_msg_l2cap_coc_send_disconnection_request_rsp_t* gecko_cmd_l2cap_coc_send_disconnection_request(uint8 connection,uint16 cid)
{
    
    gecko_cmd_msg->data.cmd_l2cap_coc_send_disconnection_request.connection=connection;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_disconnection_request.cid=cid;
    gecko_cmd_msg->header=(gecko_cmd_l2cap_coc_send_disconnection_request_id+(((3)&0xff)<<8)+(((3)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_l2cap_coc_send_disconnection_request;
}

/** 
*
* gecko_cmd_l2cap_coc_send_data
*
* Send data to a L2CAP LE connection-oriented channel. 
*
* @param connection   Handle of the LE connection for sending data
* @param cid   The CID represents the destination channel endpoint of the device sending
*  data.
* @param data_len   Array length
* @param data_data   Data to be sent. Data length must be within the range of destination channel
*  endpoint's MPS value.
*
**/

static inline struct gecko_msg_l2cap_coc_send_data_rsp_t* gecko_cmd_l2cap_coc_send_data(uint8 connection,uint16 cid,uint8 data_len, const uint8* data_data)
{
    if ((uint16_t)data_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_l2cap_coc_send_data.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_l2cap_coc_send_data;
    }

    
    gecko_cmd_msg->data.cmd_l2cap_coc_send_data.connection=connection;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_data.cid=cid;
    gecko_cmd_msg->data.cmd_l2cap_coc_send_data.data.len=data_len;
    memcpy(gecko_cmd_msg->data.cmd_l2cap_coc_send_data.data.data,data_data,data_len);
    gecko_cmd_msg->header=(gecko_cmd_l2cap_coc_send_data_id+(((4+data_len)&0xff)<<8)+(((4+data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_l2cap_coc_send_data;
}

/** 
*
* gecko_cmd_cte_transmitter_enable_connection_cte
*
* Enable different types of CTE responses on a connection. CTE response will be
* sent once requested by the peer device using the CTE Request procedure. 
*
* @param connection   Connection handle
* @param cte_types   CTE types. Bitmask of the following:
*  
*      Bit 0: AoA CTE response  
*      Bit 1: AoD CTE response with 1 us slots  
*      Bit 2: AoD CTE response with 2 us slots
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over.
*
**/

static inline struct gecko_msg_cte_transmitter_enable_connection_cte_rsp_t* gecko_cmd_cte_transmitter_enable_connection_cte(uint8 connection,uint8 cte_types,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 3)
    {
        gecko_rsp_msg->data.rsp_cte_transmitter_enable_connection_cte.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_transmitter_enable_connection_cte;
    }

    
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connection_cte.connection=connection;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connection_cte.cte_types=cte_types;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connection_cte.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_transmitter_enable_connection_cte.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_enable_connection_cte_id+(((3+switching_pattern_len)&0xff)<<8)+(((3+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_enable_connection_cte;
}

/** 
*
* gecko_cmd_cte_transmitter_disable_connection_cte
*
* Disable CTE responses on a connection. 
*
* @param connection   Connection handle
*
**/

static inline struct gecko_msg_cte_transmitter_disable_connection_cte_rsp_t* gecko_cmd_cte_transmitter_disable_connection_cte(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_cte_transmitter_disable_connection_cte.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_disable_connection_cte_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_disable_connection_cte;
}

/** 
*
* gecko_cmd_cte_transmitter_enable_connectionless_cte
*
* Start connectionless CTE transmit. CTEs will be transmitted in periodic
* advertisement packets. As a result, a periodic advertising has to be started
* prior this command. 
*
* @param handle   Periodic advertising handle
* @param cte_length   CTE length in 8 us units.
*  
*      Range: 0x02 to 0x14  
*      Time Range: 16 us to 160 us
* @param cte_type   CTE type
*  
*      0: AoA CTE  
*      1: AoD CTE with 1 us slots  
*      2: AoD CTE with 2 us slots
* @param cte_count   The number of CTEs to be transmitted in each periodic advertising interval
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over.
*
**/

static inline struct gecko_msg_cte_transmitter_enable_connectionless_cte_rsp_t* gecko_cmd_cte_transmitter_enable_connectionless_cte(uint8 handle,uint8 cte_length,uint8 cte_type,uint8 cte_count,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 5)
    {
        gecko_rsp_msg->data.rsp_cte_transmitter_enable_connectionless_cte.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_transmitter_enable_connectionless_cte;
    }

    
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connectionless_cte.handle=handle;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connectionless_cte.cte_length=cte_length;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connectionless_cte.cte_type=cte_type;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connectionless_cte.cte_count=cte_count;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_connectionless_cte.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_transmitter_enable_connectionless_cte.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_enable_connectionless_cte_id+(((5+switching_pattern_len)&0xff)<<8)+(((5+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_enable_connectionless_cte;
}

/** 
*
* gecko_cmd_cte_transmitter_disable_connectionless_cte
*
* Stop the connectionless CTE transmit. 
*
* @param handle   Periodic advertising handle
*
**/

static inline struct gecko_msg_cte_transmitter_disable_connectionless_cte_rsp_t* gecko_cmd_cte_transmitter_disable_connectionless_cte(uint8 handle)
{
    
    gecko_cmd_msg->data.cmd_cte_transmitter_disable_connectionless_cte.handle=handle;
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_disable_connectionless_cte_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_disable_connectionless_cte;
}

/** 
*
* gecko_cmd_cte_transmitter_set_dtm_parameters
*
* Set the CTE-related parameters of the LE transmitter test. 
*
* @param cte_length   Length of the Constant Tone Extension in 8 us units
*  
*      0: No CTE  
*      0x02 to 0x14: CTE length  
*  
*  Default: 0 (no CTE)
* @param cte_type   CTE type
*  
*      0: AoA CTE  
*      1: AoD CTE with 1 us slots  
*      2: AoD CTE with 2 us slots  
*  
*  Default: 0
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over. Default is the empty array.
*
**/

static inline struct gecko_msg_cte_transmitter_set_dtm_parameters_rsp_t* gecko_cmd_cte_transmitter_set_dtm_parameters(uint8 cte_length,uint8 cte_type,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 3)
    {
        gecko_rsp_msg->data.rsp_cte_transmitter_set_dtm_parameters.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_transmitter_set_dtm_parameters;
    }

    
    gecko_cmd_msg->data.cmd_cte_transmitter_set_dtm_parameters.cte_length=cte_length;
    gecko_cmd_msg->data.cmd_cte_transmitter_set_dtm_parameters.cte_type=cte_type;
    gecko_cmd_msg->data.cmd_cte_transmitter_set_dtm_parameters.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_transmitter_set_dtm_parameters.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_set_dtm_parameters_id+(((3+switching_pattern_len)&0xff)<<8)+(((3+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_set_dtm_parameters;
}

/** 
*
* gecko_cmd_cte_transmitter_clear_dtm_parameters
*
* Clear CTE-related parameters that were previously set for LE transmitter test.
* Default values will be restored for these parameters. 
*
*
**/

static inline struct gecko_msg_cte_transmitter_clear_dtm_parameters_rsp_t* gecko_cmd_cte_transmitter_clear_dtm_parameters()
{
    
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_clear_dtm_parameters_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_clear_dtm_parameters;
}

/** 
*
* gecko_cmd_cte_transmitter_enable_silabs_cte
*
* Enable Silicon Labs CTE transmit. CTEs will be transmitted in extended
* advertisement packets. As a result, extended advertising has to be started
* prior this command. 
*
* @param handle   Advertising handle
* @param cte_length   CTE length in 8 us units.
*  
*      Range: 0x02 to 0x14  
*      Time Range: 16 us to 160 us
* @param cte_type   CTE type
*  
*      0: AoA CTE  
*      1: AoD CTE with 1 us slots  
*      2: AoD CTE with 2 us slots
* @param cte_count   The number of CTEs to be transmitted in each extended advertising interval.
*  Currently only cte_count = 1 is supported.
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over.
*
**/

static inline struct gecko_msg_cte_transmitter_enable_silabs_cte_rsp_t* gecko_cmd_cte_transmitter_enable_silabs_cte(uint8 handle,uint8 cte_length,uint8 cte_type,uint8 cte_count,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 5)
    {
        gecko_rsp_msg->data.rsp_cte_transmitter_enable_silabs_cte.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_transmitter_enable_silabs_cte;
    }

    
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_silabs_cte.handle=handle;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_silabs_cte.cte_length=cte_length;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_silabs_cte.cte_type=cte_type;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_silabs_cte.cte_count=cte_count;
    gecko_cmd_msg->data.cmd_cte_transmitter_enable_silabs_cte.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_transmitter_enable_silabs_cte.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_enable_silabs_cte_id+(((5+switching_pattern_len)&0xff)<<8)+(((5+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_enable_silabs_cte;
}

/** 
*
* gecko_cmd_cte_transmitter_disable_silabs_cte
*
* Disable Silicon Labs CTE transmit. 
*
* @param handle   Advertising handle
*
**/

static inline struct gecko_msg_cte_transmitter_disable_silabs_cte_rsp_t* gecko_cmd_cte_transmitter_disable_silabs_cte(uint8 handle)
{
    
    gecko_cmd_msg->data.cmd_cte_transmitter_disable_silabs_cte.handle=handle;
    gecko_cmd_msg->header=(gecko_cmd_cte_transmitter_disable_silabs_cte_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_transmitter_disable_silabs_cte;
}

/** 
*
* gecko_cmd_cte_receiver_configure
*
* Configure the CTE sampling mode. 
*
* @param flags   Values:
*  
*      0: Disable raw sample mode, only picked IQ samples are reported (1 IQ
*      sample pair / slot)  
*      1: Enable raw sample mode, every IQ sample is reported.  
*  
*  Default: 0
*
**/

static inline struct gecko_msg_cte_receiver_configure_rsp_t* gecko_cmd_cte_receiver_configure(uint8 flags)
{
    
    gecko_cmd_msg->data.cmd_cte_receiver_configure.flags=flags;
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_configure_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_configure;
}

/** 
*
* gecko_cmd_cte_receiver_enable_connection_cte
*
* Start IQ samplings on a connection. A CTE requests will be initiated
* periodically on the given connection and IQ sampling will be made on the
* received CTE responses. 
*
* @param connection   Connection handle
* @param interval   Measurement interval
*  
*      0: No interval. The request is initiated only once.  
*      Other values N: Initiate the request every N-th connection events
* @param cte_length   Minimum CTE length requested in 8 us units.
*  
*      Range: 0x02 to 0x14  
*      Time Range: 16 us to 160 us
* @param cte_type   Requested CTE type
*  
*      0: AoA CTE  
*      1: AoD CTE with 1 us slots  
*      2: AoD CTE with 2 us slots
* @param slot_durations   Slot durations
*  
*      1: Switching and sampling slots are 1 us each  
*      2: Switching and sampling slots are 2 us each
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over.
*
* Events generated
*
* gecko_evt_cte_receiver_connection_iq_report - Triggered when IQ samples have been received.
*
**/

static inline struct gecko_msg_cte_receiver_enable_connection_cte_rsp_t* gecko_cmd_cte_receiver_enable_connection_cte(uint8 connection,uint16 interval,uint8 cte_length,uint8 cte_type,uint8 slot_durations,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 7)
    {
        gecko_rsp_msg->data.rsp_cte_receiver_enable_connection_cte.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_receiver_enable_connection_cte;
    }

    
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.connection=connection;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.interval=interval;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.cte_length=cte_length;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.cte_type=cte_type;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.slot_durations=slot_durations;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_receiver_enable_connection_cte.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_enable_connection_cte_id+(((7+switching_pattern_len)&0xff)<<8)+(((7+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_enable_connection_cte;
}

/** 
*
* gecko_cmd_cte_receiver_disable_connection_cte
*
* Stop the IQ sampling on a connection. CTEs will not be requested on the given
* connection. 
*
* @param connection   Connection handle
*
**/

static inline struct gecko_msg_cte_receiver_disable_connection_cte_rsp_t* gecko_cmd_cte_receiver_disable_connection_cte(uint8 connection)
{
    
    gecko_cmd_msg->data.cmd_cte_receiver_disable_connection_cte.connection=connection;
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_disable_connection_cte_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_disable_connection_cte;
}

/** 
*
* gecko_cmd_cte_receiver_enable_connectionless_cte
*
* Start IQ sampling on a periodic advertising synchronization. IQ samples are
* taken on each CTE found in the periodic advertisements. 
*
* @param sync   Periodic advertising synchronization handle
* @param slot_durations   Slot durations
*  
*      1: Switching and sampling slots are 1 us each  
*      2: Switching and sampling slots are 2 us each
* @param cte_count   * 0: Sample and report all available CTEs  
*      Other values: Maximum number of sampled CTEs in each periodic advertising
*      interval
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over.
*
* Events generated
*
* gecko_evt_cte_receiver_connectionless_iq_report - Triggered when IQ samples have been received.
*
**/

static inline struct gecko_msg_cte_receiver_enable_connectionless_cte_rsp_t* gecko_cmd_cte_receiver_enable_connectionless_cte(uint8 sync,uint8 slot_durations,uint8 cte_count,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_cte_receiver_enable_connectionless_cte.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_receiver_enable_connectionless_cte;
    }

    
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connectionless_cte.sync=sync;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connectionless_cte.slot_durations=slot_durations;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connectionless_cte.cte_count=cte_count;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_connectionless_cte.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_receiver_enable_connectionless_cte.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_enable_connectionless_cte_id+(((4+switching_pattern_len)&0xff)<<8)+(((4+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_enable_connectionless_cte;
}

/** 
*
* gecko_cmd_cte_receiver_disable_connectionless_cte
*
* Stop IQ sampling on a periodic advertising synchronization. 
*
* @param sync   Periodic advertising synchronization handle
*
**/

static inline struct gecko_msg_cte_receiver_disable_connectionless_cte_rsp_t* gecko_cmd_cte_receiver_disable_connectionless_cte(uint8 sync)
{
    
    gecko_cmd_msg->data.cmd_cte_receiver_disable_connectionless_cte.sync=sync;
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_disable_connectionless_cte_id+(((1)&0xff)<<8)+(((1)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_disable_connectionless_cte;
}

/** 
*
* gecko_cmd_cte_receiver_set_dtm_parameters
*
* Set CTE-related parameters of LE receiver test. 
*
* @param cte_length   Expected CTE length in 8 us units
*  
*      0: No CTE  
*      0x02 to 0x14: Expected CTE length  
*  
*  Default: 0 (no CTE)
* @param cte_type   Expected CTE type
*  
*      0: Expect AoA CTE  
*      1: Expect AoD CTE with 1 us slots  
*      2: Expect AoD CTE with 2 us slots  
*  
*  Default: 0
* @param slot_durations   Slot durations
*  
*      1: Switching and sampling slots are 1 us each  
*      2: Switching and sampling slots are 2 us each  
*  
*  Default: 1
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over. Default: empty array
*
* Events generated
*
* gecko_evt_cte_receiver_dtm_iq_report - Triggered when IQ samples have been received.
*
**/

static inline struct gecko_msg_cte_receiver_set_dtm_parameters_rsp_t* gecko_cmd_cte_receiver_set_dtm_parameters(uint8 cte_length,uint8 cte_type,uint8 slot_durations,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 4)
    {
        gecko_rsp_msg->data.rsp_cte_receiver_set_dtm_parameters.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_receiver_set_dtm_parameters;
    }

    
    gecko_cmd_msg->data.cmd_cte_receiver_set_dtm_parameters.cte_length=cte_length;
    gecko_cmd_msg->data.cmd_cte_receiver_set_dtm_parameters.cte_type=cte_type;
    gecko_cmd_msg->data.cmd_cte_receiver_set_dtm_parameters.slot_durations=slot_durations;
    gecko_cmd_msg->data.cmd_cte_receiver_set_dtm_parameters.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_receiver_set_dtm_parameters.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_set_dtm_parameters_id+(((4+switching_pattern_len)&0xff)<<8)+(((4+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_set_dtm_parameters;
}

/** 
*
* gecko_cmd_cte_receiver_clear_dtm_parameters
*
* Clear CTE-related parameters that were previously set for LE receiver test.
* Default values will be restored for these parameters. 
*
*
**/

static inline struct gecko_msg_cte_receiver_clear_dtm_parameters_rsp_t* gecko_cmd_cte_receiver_clear_dtm_parameters()
{
    
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_clear_dtm_parameters_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_clear_dtm_parameters;
}

/** 
*
* gecko_cmd_cte_receiver_enable_silabs_cte
*
* Enable IQ sampling of Silicon Labs CTE found in extended advertisements. 
*
* @param slot_durations   Slot durations
*  
*      1: Switching and sampling slots are 1 us each  
*      2: Switching and sampling slots are 2 us each
* @param cte_count   * 0: Sample and report all available CTEs  
*      Other values: Maximum number of sampled CTEs in each extended advertising
*      interval
* @param switching_pattern_len   Array length
* @param switching_pattern_data   Antenna switching pattern. Antennas will be switched in this order with the
*  antenna switch pins during CTE. If the CTE is longer than the switching
*  pattern, the pattern starts over.
*
* Events generated
*
* gecko_evt_cte_receiver_silabs_iq_report - Triggered when IQ samples of Silicon Labs CTE have been received.
*
**/

static inline struct gecko_msg_cte_receiver_enable_silabs_cte_rsp_t* gecko_cmd_cte_receiver_enable_silabs_cte(uint8 slot_durations,uint8 cte_count,uint8 switching_pattern_len, const uint8* switching_pattern_data)
{
    if ((uint16_t)switching_pattern_len > BGLIB_MSG_MAX_PAYLOAD - 3)
    {
        gecko_rsp_msg->data.rsp_cte_receiver_enable_silabs_cte.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_cte_receiver_enable_silabs_cte;
    }

    
    gecko_cmd_msg->data.cmd_cte_receiver_enable_silabs_cte.slot_durations=slot_durations;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_silabs_cte.cte_count=cte_count;
    gecko_cmd_msg->data.cmd_cte_receiver_enable_silabs_cte.switching_pattern.len=switching_pattern_len;
    memcpy(gecko_cmd_msg->data.cmd_cte_receiver_enable_silabs_cte.switching_pattern.data,switching_pattern_data,switching_pattern_len);
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_enable_silabs_cte_id+(((3+switching_pattern_len)&0xff)<<8)+(((3+switching_pattern_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_enable_silabs_cte;
}

/** 
*
* gecko_cmd_cte_receiver_disable_silabs_cte
*
* Disable IQ sampling of Silicon Labs CTE. 
*
*
**/

static inline struct gecko_msg_cte_receiver_disable_silabs_cte_rsp_t* gecko_cmd_cte_receiver_disable_silabs_cte()
{
    
    gecko_cmd_msg->header=(gecko_cmd_cte_receiver_disable_silabs_cte_id+(((0)&0xff)<<8)+(((0)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_cte_receiver_disable_silabs_cte;
}

/** 
*
* gecko_cmd_user_message_to_target
*
* Used by an NCP host to send a message to the target application on device. The
* application on target is must send the response with
* gecko_send_rsp_user_message_to_target. 
*
* @param data_len   Array length
* @param data_data   The message
*
**/

static inline struct gecko_msg_user_message_to_target_rsp_t* gecko_cmd_user_message_to_target(uint8 data_len, const uint8* data_data)
{
    if ((uint16_t)data_len > BGLIB_MSG_MAX_PAYLOAD - 1)
    {
        gecko_rsp_msg->data.rsp_user_message_to_target.result = bg_err_command_too_long;
        return &gecko_rsp_msg->data.rsp_user_message_to_target;
    }

    
    gecko_cmd_msg->data.cmd_user_message_to_target.data.len=data_len;
    memcpy(gecko_cmd_msg->data.cmd_user_message_to_target.data.data,data_data,data_len);
    gecko_cmd_msg->header=(gecko_cmd_user_message_to_target_id+(((1+data_len)&0xff)<<8)+(((1+data_len)&0x700)>>8));
    
    gecko_handle_command(gecko_cmd_msg->header,&gecko_cmd_msg->data.payload);
    
    return &gecko_rsp_msg->data.rsp_user_message_to_target;
}
#ifdef __cplusplus
}
#endif

#endif
