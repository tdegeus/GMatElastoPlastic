/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlastic

*/

#ifndef GMATELASTOPLASTIC_CONFIG_H
#define GMATELASTOPLASTIC_CONFIG_H

#ifdef GMATELASTOPLASTIC_ENABLE_ASSERT

    #define GMATELASTOPLASTIC_ASSERT(expr) GMATELASTOPLASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)
    #define GMATELASTOPLASTIC_ASSERT_IMPL(expr, file, line) \
        if (!(expr)) { \
            throw std::runtime_error( \
                std::string(file) + ':' + std::to_string(line) + \
                ": assertion failed (" #expr ") \n\t"); \
        }

#else

    #define GMATELASTOPLASTIC_ASSERT(expr)

#endif

#endif
