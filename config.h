/*
 Copyright (c)  2016 Afanasieva A.
                     Bekbulatov R.
                     Nazarov V.
                     Halikov P.
 Last source on https://github.com/uburuntu/diff_equations
*/

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
 */

/*
*
*     Area 10 exercise:
*      _____________________
*    ->          |          |
*     |          |          |
*    ->    02    |    12    |
*     |          |          |
*    ->__ __ __ _|__________|
*    ->          |          |
*     |          |          |
*    ->    01    |    11    |
*     |          |          |
*    ->__ __ __ _|__________|_ __ __ __
*                |          |         ->
*                |          |         ->
*                |    10    |    20   ->
*                |          |         ->
*                |__________|_________->
*
*/

/*
*
*     Area 11 exercise:
*      ________________________________
*     |          |          |          |
*     |                                |
*     |    01    |    11    |    21    |
*     |                                |
*     |__ __ __ _|__________|_ __ __ __|
*     |          |          |          |
*     |          |          |          |
*     |    00    |    10    |    20    |
*     |          |          |          |
*     |__________|          |__________|
*
*/
