/// First prog config

#define RHO_0   0.5
#define RHO_G   0.3
#define W       0.1

#define STAT_SOL_EPS    1.e-12

typedef enum enum_calc_type_t
{
  SMOOTH,
  NO_SMOOTH,
  STAT_SOL_SRCH,
  EIG_FUNC_INIT
} calc_type_t;

// Change this for different calculation type
static const calc_type_t calc_type = STAT_SOL_SRCH;



/// Common config (for first program and eigen)

#define A_LENGHT 1.
#define C_LENGHT 2.
#define B_LENGHT 3.

typedef enum enum_grid_type_t
{
  SQUARE,
  VOLODYA_9,
  RAMZAN_10,
  NASTYA_11
} grid_type_t;

// Change this for different grids
static const grid_type_t grid_type = VOLODYA_9;

/*
 *
 *     Area 9 exercise:
 *     ________________________________
 *     |         |          |          |
 *     |         |          |          |
 *     |   02    |    12    |    22    |
 *     |         |          |          |
 *     |_________|__________|__________|
 *     |         |          |          |
 *     |         |          |          |
 *     |   01    |    11    |    21    |
 *     |         |          |          |
 *     |_________|__________|__________|
 *                          |          |
 *                          |          |
 *         00         10    |    20    |
 *                          |          |
 *                          |__________|
 *
 *     <-------------------><---------->
 *               c               a
 *     <------------------------------->
 *                     b
 *
 *
 *  Statuses:
 *  0 - inner                                  +
 *  1 - x = 0, y = (a, b) && x = c, y = (0, a) modified
 *  2 - x = b, y = (0, b)                      +
 *  3 - y = 0, x = (c, b) && y = a, x = (0, c) modified
 *  4 - y = b, x = (0, b)                      +
 *  5 - x = 0, y = a && x = c, y = 0           modified
 *  6 - x = b, y = 0                           +
 *  7 - x = 0. y = b                           +
 *  8 - x = b, y = b                           +
 *
 */
