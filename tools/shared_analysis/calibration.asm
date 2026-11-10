; For n >= 0 this executes exactly 3*n + 4 instructions and returns n.
; Used only to verify the optional instruction-counting client.
.code
icount_loop PROC
    xor eax, eax
    test ecx, ecx
    jz finished
again:
    inc eax
    dec ecx
    jnz again
finished:
    ret
icount_loop ENDP

; An application-owned no-op boundary, intercepted by the optional client.
; Arguments: operation (1 = begin, 0 = end), NUL-terminated region label.
icount_marker PROC
    ret
icount_marker ENDP
END
