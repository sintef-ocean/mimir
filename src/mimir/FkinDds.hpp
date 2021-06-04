#pragma once
/**
   @file FkinDds.hpp
   @brief Convenience header for FKIN DDS interface and DDS header files.
*/
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif
#include "FKIN/fkin_types_DCPS.hpp"
#include "Ratatosk/basic_types_DCPS.hpp"
#include <dds/domain/DomainParticipant.hpp>
#include <dds/pub/Publisher.hpp>
#include <dds/sub/Subscriber.hpp>
#ifdef _MSC_VER
#pragma warning(pop)
#endif
