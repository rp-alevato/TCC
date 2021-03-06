#ifndef RETARGETSERIALHALCONFIG_H
#define RETARGETSERIALHALCONFIG_H

#include "hal-config.h"

#if BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_USART0
// USART0
#define RETARGET_UART       USART0
#define RETARGET_CLK        cmuClock_USART0
#define RETARGET_UART_INDEX 0
#if defined(_SILICON_LABS_GECKO_INTERNAL_SDID_103)
// EFM32TG11 has merged USART interrupts
#define RETARGET_IRQ_NAME   USART0_IRQHandler
#define RETARGET_IRQn       USART0_IRQn
#else
#define RETARGET_IRQ_NAME   USART0_RX_IRQHandler
#define RETARGET_IRQn       USART0_RX_IRQn
#endif
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_USART1
// USART1
#define RETARGET_UART       USART1
#define RETARGET_CLK        cmuClock_USART1
#define RETARGET_UART_INDEX 1
#if defined(_SILICON_LABS_GECKO_INTERNAL_SDID_103)
// EFM32TG11 has merged USART interrupts
#define RETARGET_IRQ_NAME   USART1_IRQHandler
#define RETARGET_IRQn       USART1_IRQn
#else
#define RETARGET_IRQ_NAME   USART1_RX_IRQHandler
#define RETARGET_IRQn       USART1_RX_IRQn
#endif
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_USART2
// USART2
#define RETARGET_UART       USART2
#define RETARGET_CLK        cmuClock_USART2
#define RETARGET_UART_INDEX 2
#if defined(_SILICON_LABS_GECKO_INTERNAL_SDID_103)
// EFM32TG11 has merged USART interrupts
#define RETARGET_IRQ_NAME   USART2_IRQHandler
#define RETARGET_IRQn       USART2_IRQn
#else
#define RETARGET_IRQ_NAME   USART2_RX_IRQHandler
#define RETARGET_IRQn       USART2_RX_IRQn
#endif
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_USART3
// USART3
#define RETARGET_UART       USART3
#define RETARGET_CLK        cmuClock_USART3
#define RETARGET_UART_INDEX 3
#if defined(_SILICON_LABS_GECKO_INTERNAL_SDID_103)
// EFM32TG11 has merged USART interrupts
#define RETARGET_IRQ_NAME   USART3_IRQHandler
#define RETARGET_IRQn       USART3_IRQn
#else
#define RETARGET_IRQ_NAME   USART3_RX_IRQHandler
#define RETARGET_IRQn       USART3_RX_IRQn
#endif
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_USART4
// USART4
#define RETARGET_UART       USART4
#define RETARGET_CLK        cmuClock_USART4
#define RETARGET_UART_INDEX 4
#define RETARGET_IRQ_NAME   USART4_RX_IRQHandler
#define RETARGET_IRQn       USART4_RX_IRQn
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_USART5
// USART5
#define RETARGET_UART       USART5
#define RETARGET_CLK        cmuClock_USART5
#define RETARGET_UART_INDEX 5
#define RETARGET_IRQ_NAME   USART5_RX_IRQHandler
#define RETARGET_IRQn       USART5_RX_IRQn
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_UART0
// UART0
#define RETARGET_UART       UART0
#define RETARGET_CLK        cmuClock_UART0
#define RETARGET_UART_INDEX 0
#if defined(_SILICON_LABS_GECKO_INTERNAL_SDID_103)
// EFM32TG11 has merged UART interrupts
#define RETARGET_IRQ_NAME   UART0_IRQHandler
#define RETARGET_IRQn       UART0_IRQn
#else
#define RETARGET_IRQ_NAME   UART0_RX_IRQHandler
#define RETARGET_IRQn       UART0_RX_IRQn
#endif
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_UART1
// UART1
#define RETARGET_UART       UART1
#define RETARGET_CLK        cmuClock_UART1
#define RETARGET_UART_INDEX 1
#if defined(_SILICON_LABS_GECKO_INTERNAL_SDID_103)
// EFM32TG11 has merged UART interrupts
#define RETARGET_IRQ_NAME   UART1_IRQHandler
#define RETARGET_IRQn       UART1_IRQn
#else
#define RETARGET_IRQ_NAME   UART1_RX_IRQHandler
#define RETARGET_IRQn       UART1_RX_IRQn
#endif
#define RETARGET_USART      1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_LEUART0
// LEUART0
#define RETARGET_UART       LEUART0
#define RETARGET_CLK        cmuClock_LEUART0
#define RETARGET_UART_INDEX 0
#define RETARGET_IRQ_NAME   LEUART0_IRQHandler
#define RETARGET_IRQn       LEUART0_IRQn
#define RETARGET_LEUART     1
#elif BSP_SERIAL_APP_PORT == HAL_SERIAL_PORT_LEUART1
// LEUART1
#define RETARGET_UART       LEUART1
#define RETARGET_CLK        cmuClock_LEUART1
#define RETARGET_UART_INDEX 1
#define RETARGET_IRQ_NAME   LEUART1_IRQHandler
#define RETARGET_IRQn       LEUART1_IRQn
#define RETARGET_LEUART     1
#endif

#if defined(RETARGET_USART)
// UART is an USART
#define RETARGET_TX         USART_Tx
#define RETARGET_RX         USART_Rx
#elif defined(RETARGET_LEUART)
// UART is a LEUART
#define RETARGET_TX         LEUART_Tx
#define RETARGET_RX         LEUART_Rx
#else
#error "Invalid UART type"
#endif

#define RETARGET_TXPORT           BSP_SERIAL_APP_TX_PORT
#define RETARGET_TXPIN            BSP_SERIAL_APP_TX_PIN
#define RETARGET_RXPORT           BSP_SERIAL_APP_RX_PORT
#define RETARGET_RXPIN            BSP_SERIAL_APP_RX_PIN

#if defined(_SILICON_LABS_32B_SERIES_0)
// Series 0 devices only have one location
#define RETARGET_LOCATION         BSP_SERIAL_APP_ROUTE_LOC
#else
// Series 1 devices have one location per pin
#define RETARGET_TX_LOCATION      BSP_SERIAL_APP_TX_LOC
#define RETARGET_RX_LOCATION      BSP_SERIAL_APP_RX_LOC
#endif

#if defined(BSP_SERIAL_APP_CTS_PORT)
// Only Series 1+ USARTs have CTS
#define RETARGET_CTSPORT          BSP_SERIAL_APP_CTS_PORT
#define RETARGET_CTSPIN           BSP_SERIAL_APP_CTS_PIN
#define RETARGET_CTS_LOCATION     BSP_SERIAL_APP_CTS_LOC
#endif
#if defined(BSP_SERIAL_APP_RTS_PORT)
// Only Series 1+ USARTs have RTS
#define RETARGET_RTSPORT          BSP_SERIAL_APP_RTS_PORT
#define RETARGET_RTSPIN           BSP_SERIAL_APP_RTS_PIN
#define RETARGET_RTS_LOCATION     BSP_SERIAL_APP_RTS_LOC
#endif

#if HAL_VCOM_ENABLE && defined(BSP_VCOM_ENABLE_PORT)
#define RETARGET_PERIPHERAL_ENABLE()    \
  GPIO_PinModeSet(BSP_VCOM_ENABLE_PORT, \
                  BSP_VCOM_ENABLE_PIN,  \
                  gpioModePushPull,     \
                  1)
#else
#define RETARGET_PERIPHERAL_ENABLE()
#endif

#endif // RETARGETSERIALHALCONFIG_H
