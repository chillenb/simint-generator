#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Initializes the OS3C2E functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_os3c2e_init(void);


/*! \brief Initializes the OS3C2E 1st derivative functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_os3c2e_deriv1_init(void);


/*! \brief Finalizes the OS3C2E functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_os3c2e_finalize(void);


/*! \brief Finalizes the OS3C2E 1st derivative functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_os3c2e_deriv1_finalize(void);


#ifdef __cplusplus
}
#endif

