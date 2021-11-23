	.file	"exchange.c"
	.text
.Ltext0:
	.globl	exchange
	.type	exchange, @function
exchange:
.LVL0:
.LFB0:
	.file 1 "exchange.c"
	.loc 1 2 1 view -0
	.cfi_startproc
	.loc 1 2 1 is_stmt 0 view .LVU1
	endbr64
	.loc 1 3 2 is_stmt 1 view .LVU2
	.loc 1 3 7 is_stmt 0 view .LVU3
	movq	(%rdi), %rax
.LVL1:
	.loc 1 4 2 is_stmt 1 view .LVU4
	.loc 1 4 6 is_stmt 0 view .LVU5
	movq	%rsi, (%rdi)
	.loc 1 5 2 is_stmt 1 view .LVU6
	.loc 1 6 1 is_stmt 0 view .LVU7
	ret
	.cfi_endproc
.LFE0:
	.size	exchange, .-exchange
	.globl	main
	.type	main, @function
main:
.LFB1:
	.loc 1 9 1 is_stmt 1 view -0
	.cfi_startproc
	endbr64
	.loc 1 11 1 view .LVU9
	ret
	.cfi_endproc
.LFE1:
	.size	main, .-main
.Letext0:
	.section	.debug_info,"",@progbits
.Ldebug_info0:
	.long	0x9f
	.value	0x4
	.long	.Ldebug_abbrev0
	.byte	0x8
	.uleb128 0x1
	.long	.LASF0
	.byte	0xc
	.long	.LASF1
	.long	.LASF2
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.long	.Ldebug_line0
	.uleb128 0x2
	.long	.LASF3
	.byte	0x1
	.byte	0x8
	.byte	0x6
	.quad	.LFB1
	.quad	.LFE1-.LFB1
	.uleb128 0x1
	.byte	0x9c
	.uleb128 0x3
	.long	.LASF4
	.byte	0x1
	.byte	0x1
	.byte	0x6
	.long	0x95
	.quad	.LFB0
	.quad	.LFE0-.LFB0
	.uleb128 0x1
	.byte	0x9c
	.long	0x95
	.uleb128 0x4
	.string	"xp"
	.byte	0x1
	.byte	0x1
	.byte	0x15
	.long	0x9c
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x4
	.string	"y"
	.byte	0x1
	.byte	0x1
	.byte	0x1f
	.long	0x9c
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x5
	.string	"x"
	.byte	0x1
	.byte	0x3
	.byte	0x7
	.long	0x95
	.long	.LLST0
	.long	.LVUS0
	.byte	0
	.uleb128 0x6
	.byte	0x8
	.byte	0x5
	.long	.LASF5
	.uleb128 0x7
	.byte	0x8
	.long	0x95
	.byte	0
	.section	.debug_abbrev,"",@progbits
.Ldebug_abbrev0:
	.uleb128 0x1
	.uleb128 0x11
	.byte	0x1
	.uleb128 0x25
	.uleb128 0xe
	.uleb128 0x13
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x1b
	.uleb128 0xe
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x10
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x2
	.uleb128 0x2e
	.byte	0
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x2117
	.uleb128 0x19
	.byte	0
	.byte	0
	.uleb128 0x3
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x2117
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x4
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x5
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x6
	.uleb128 0x24
	.byte	0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.byte	0
	.byte	0
	.uleb128 0x7
	.uleb128 0xf
	.byte	0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.byte	0
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
.LVUS0:
	.uleb128 .LVU4
	.uleb128 0
.LLST0:
	.quad	.LVL1-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x1
	.byte	0x50
	.quad	0
	.quad	0
	.section	.debug_aranges,"",@progbits
	.long	0x2c
	.value	0x2
	.long	.Ldebug_info0
	.byte	0x8
	.byte	0
	.value	0
	.value	0
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.quad	0
	.quad	0
	.section	.debug_line,"",@progbits
.Ldebug_line0:
	.section	.debug_str,"MS",@progbits,1
.LASF4:
	.string	"exchange"
.LASF2:
	.string	"/home/sunshine/proj/c/test_disassemble"
.LASF1:
	.string	"exchange.c"
.LASF5:
	.string	"long int"
.LASF0:
	.string	"GNU C17 9.3.0 -mtune=generic -march=x86-64 -g -Og -fasynchronous-unwind-tables -fstack-protector-strong -fstack-clash-protection -fcf-protection"
.LASF3:
	.string	"main"
	.ident	"GCC: (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
