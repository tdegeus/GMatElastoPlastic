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

#define GMATELASTOPLASTIC_VERSION_MAJOR 0
#define GMATELASTOPLASTIC_VERSION_MINOR 2
#define GMATELASTOPLASTIC_VERSION_PATCH 0

#define GMATELASTOPLASTIC_VERSION_AT_LEAST(x,y,z) \
    (GMATELASTOPLASTIC_VERSION_MAJOR > x || (GMATELASTOPLASTIC_VERSION_MAJOR >= x && \
    (GMATELASTOPLASTIC_VERSION_MINOR > y || (GMATELASTOPLASTIC_VERSION_MINOR >= y && \
                                             GMATELASTOPLASTIC_VERSION_PATCH >= z))))

#endif
