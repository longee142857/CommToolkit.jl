## ALINX AX7020（Zynq-7020 / XC7Z020-CLG400）PL 占位约束
##
## 引脚摘自黑金 AX7020 用户手册里常见的 PL 分配，不是本仓库实测的结果。
## 不同批次、核心板和底板请对照自己的原理图再改。本目录没有上板记录。
##
##   PL_GCLK  50 MHz  单端  U18     这是 PL 晶振，不是 PS 的 33.333 MHz（E7）
##   KEY1     低有效  N15     复位
##   KEY2     低有效  N16     切换 Eb/N0
##   LED1     低有效  M14     led[0] 心跳
##   LED2     低有效  M15     led[1] snr_index[0]
##   LED3     低有效  K16     led[2] snr_index[1]
##   LED4     低有效  J16     led[3] snr_index[2]
##
## BANK 电压按手册写 LVCMOS33。若某只脚所在 BANK 不是 3.3 V，改 IOSTANDARD，
## 不要照抄。

set_property PACKAGE_PIN U18 [get_ports clk]
set_property IOSTANDARD LVCMOS33 [get_ports clk]
create_clock -period 20.000 -name pl_clk [get_ports clk]

set_property PACKAGE_PIN N15 [get_ports rst_n]
set_property IOSTANDARD LVCMOS33 [get_ports rst_n]

set_property PACKAGE_PIN N16 [get_ports key_n]
set_property IOSTANDARD LVCMOS33 [get_ports key_n]

set_property PACKAGE_PIN M14 [get_ports {led[0]}]
set_property PACKAGE_PIN M15 [get_ports {led[1]}]
set_property PACKAGE_PIN K16 [get_ports {led[2]}]
set_property PACKAGE_PIN J16 [get_ports {led[3]}]
set_property IOSTANDARD LVCMOS33 [get_ports {led[*]}]
