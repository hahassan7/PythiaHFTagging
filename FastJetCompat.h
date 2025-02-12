#pragma once

// Compatibility layer for FastJet's auto_ptr issue
#ifdef __cplusplus
#if __cplusplus >= 201402L
#define FASTJET_COMPAT
#include <memory>
namespace fastjet {
    namespace internal {
        template<typename T>
        using auto_ptr = std::unique_ptr<T>;
    }
}
namespace std {
    template<typename T>
    using auto_ptr = std::unique_ptr<T>;
}
#endif
#endif 