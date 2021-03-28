#if 0
int usleep(useconds_t) __asm("_" "usleep" );
#endif
#define USLEEP_RET_TYPE int
#define USLEEP_RET_VAL 0
#define USLEEP_ARG_TYPE useconds_t
