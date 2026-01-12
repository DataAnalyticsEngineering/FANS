//
// Created by Alex Hocks on 09/01/26.
//

#ifndef SERIALIZATION_H
#define SERIALIZATION_H
#include "logging.h"

#include <cstdint>
#include <vector>
#include <type_traits>

template<typename T, typename Base = std::vector<T>, typename r = std::conditional_t<sizeof(T)<=4, std::true_type, void>>
class SerializationBuffer : public Base {
    // Buffer always has size multiple of 4 byte
public:
    using Base::Base;
    using size_type = typename Base::size_type;
    void resize(size_type new_size)
    {
        const auto num_bytes = new_size * data_size;
        const auto delta_bytes = (4 - (num_bytes % 4)) % 4;
        new_size += delta_bytes / data_size;
        Base::resize(new_size);
    }

    void resize(size_type new_size, const T& value)
    {
        const auto num_bytes = new_size * data_size;
        const auto delta_bytes = (4 - (num_bytes % 4)) % 4;
        new_size += delta_bytes / data_size;
        Base::resize(new_size, value);
    }

    void add_size(size_type delta) { resize(Base::size() + delta); }

    std::vector<int>& as_int_vector() { return *reinterpret_cast<std::vector<int>*>(this); }

    static SerializationBuffer& from_int_vector(std::vector<int>& vec) { return *reinterpret_cast<SerializationBuffer*>(&vec); }

private:
    static constexpr auto data_size = sizeof(T);
};

class Serializable {
    /**
     * Serializes any object that implements this.
     * Buffer elements will be bytes, but total count will always be padded to be multiple of 4.
     * Data is run length encoded.
     */
public:
    // run length type
    using length_t = uint64_t;
    // run length block size
    static constexpr auto header_len = sizeof(length_t);
    // byte data type
    using data_t   = uint8_t;
    // buffer type
    using buffer_t = SerializationBuffer<data_t>;
    // registry type: pointer to data, data size, true: run length encoded, false: direct read write
    using reg_entry_t = std::tuple<void*, length_t, bool>;
    using registry_t = std::vector<reg_entry_t>;

    virtual ~Serializable() = default;

    /**
     * Serializes this object
     * @return serial buffer
     */
    buffer_t serialize_full()
    {
        if (manual_deserialize) {
            buffer_t buffer;
            serialize_override(buffer, 0);
            return buffer;
        }

        Log::io->trace() << "Serialization::serialize_full() - register_serialization\n";
        registry_t registry;
        register_serialization(registry);

        Log::io->trace() << "Serialization::serialize_full() - accumulate\n";
        const length_t total_size = std::accumulate(registry.begin(), registry.end(), 0UL,
            [](const length_t acc, const reg_entry_t &e) {
            const bool use_run_length = std::get<2>(e);
            length_t size = std::get<1>(e);
            if (use_run_length) size += header_len;
            return acc + size;
        });

        buffer_t buffer;
        buffer.resize(total_size);
        Log::io->trace() << "Serialization::serialize_full() - buffer size=" << total_size << " addr=" << static_cast<void*>(std::data(buffer)) << "\n";

        Log::io->trace() << "Serialization::serialize_full() - process entries\n";
        length_t offset = 0UL;
        length_t pos = 0UL;
        for (const reg_entry_t &e : registry) {
            const bool use_run_length = std::get<2>(e);
            length_t size             = std::get<1>(e);
            void* data                = std::get<0>(e);
            Log::io->trace() << "Serialization::serialize_full() - entry: pos=" << pos << " addr=" << data << " offset=" << offset << " useRL=" << use_run_length << " size=" << size << "\n";

            if (use_run_length) {
                *reinterpret_cast<length_t*>(std::data(buffer) + offset) = size;
                offset += header_len;
            }
            std::memcpy(std::data(buffer) + offset, data, size);
            offset += size;
            pos++;
        }

        Log::io->trace() << "Serialization::serialize_full() - done\n";
        return buffer;
    }

    virtual void register_serialization(registry_t &r) {}

    virtual void init_deserialization() {}

    virtual registry_t::iterator late_init_deserialization(registry_t &r, registry_t::iterator begin) { return r.begin(); }

    bool manual_deserialize = false;
    bool manual_serialize = false;

    /**
     * Deserialize this object using provided data
     * @param buffer
     */
    void deserialize_full(buffer_t &buffer)
    {
        if (manual_deserialize) {
            deserialize_override(buffer, 0);
            return;
        }

        init_deserialization(); // alloc mem if needed

        // gather serialization description
        registry_t registry;
        register_serialization(registry);
        // place correct lengths into descriptor list, etc
        late_init_deserialization(registry, registry.begin());

        length_t offset = 0UL;
        for (reg_entry_t &e : registry) {
            const bool use_run_length = std::get<2>(e);
            length_t &size            = std::get<1>(e);
            void* data                = std::get<0>(e);

            if (use_run_length) {
                size = *reinterpret_cast<const length_t*>(std::data(buffer) + offset);
                offset += header_len;
            }
            std::memcpy(data, std::data(buffer) + offset, size);
            offset += size;
        }
    }

    virtual length_t deserialize_override(buffer_t &buffer, length_t offset) { return offset; }
    virtual length_t serialize_override(buffer_t &buffer, length_t offset) { return offset; }
};

#endif
