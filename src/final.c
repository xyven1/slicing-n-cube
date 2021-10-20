#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

long get_file_size(const char *path) {
  struct stat stat_buf;
  int rc = stat(path, &stat_buf);
  return (rc == 0) ? stat_buf.st_size : -1;
}

size_t read_from_file(const char *path, char **buf_ptr) {
  const long maybe_file_size = get_file_size(path);
  if (maybe_file_size <= 0) {
    return 0;
  }
  const size_t file_size = (unsigned long)maybe_file_size;
  char *buf = (char *)malloc(file_size);
  if (buf == NULL) {
    return 0;
  }
  FILE *f = fopen(path, "rb");
  if (f == NULL) {
    free(buf);
    return 0;
  }
  const size_t elements_read = fread(buf, sizeof(char), file_size, f);
  const int close = fclose(f);
  if (elements_read != file_size || close != 0) {
    free(buf);
    return 0;
  }
  *buf_ptr = buf;
  return file_size;
}

void printb(const char *buf, size_t num_bytes) {
  for (size_t i = 0; i < num_bytes; ++i) {
    char str[9];
    str[0] = (buf[i] & 0x80) ? '1' : '0';
    str[1] = (buf[i] & 0x40) ? '1' : '0';
    str[2] = (buf[i] & 0x20) ? '1' : '0';
    str[3] = (buf[i] & 0x10) ? '1' : '0';
    str[4] = (buf[i] & 0x08) ? '1' : '0';
    str[5] = (buf[i] & 0x04) ? '1' : '0';
    str[6] = (buf[i] & 0x02) ? '1' : '0';
    str[7] = (buf[i] & 0x01) ? '1' : '0';
    str[8] = '\0';
    printf("%s", str);
  }
}

int get_leading_zeros(const char *buf, size_t num_bytes) {
  for (size_t i = 0; i < num_bytes; ++i) {
    for (int j = 7; j >= 0; --j) {
      if (buf[i] & (1 << j)) {
        return (int)(i * 8) + 7 - j;
      }
    }
  }
  return (int)(num_bytes * 8);
}

int get_leading_ones(const char *buf, size_t num_bytes) {
  for (size_t i = 0; i < num_bytes; ++i) {
    for (int j = 7; j >= 0; --j) {
      if (!(buf[i] & (1 << j))) {
        return (int)(i * 8) + 7 - j;
      }
    }
  }
  return (int)(num_bytes * 8);
}

int combine_usr_mss(const char *usr, size_t usr_len, const char *mss,
                    size_t mss_len) {
  for (size_t i = 0; i < usr_len; i += 10) {
    const int leading_zeros = get_leading_zeros(usr + i, 10);
    const uint64_t mask =
        (leading_zeros < 64) ? 0xFFFFFFFFFFFFFFFF >> leading_zeros : 0;
    for (size_t j = mss_len - 10; j < mss_len; j -= 10) {
      const uint64_t usr_a = *(const uint64_t *)(usr + i);
      const uint64_t mss_a = *(const uint64_t *)(mss + j);
      const uint16_t usr_b = *(const uint16_t *)(usr + i + 8);
      const uint16_t mss_b = *(const uint16_t *)(mss + j + 8);
      #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        if ((mss_a | __builtin_bswap64(mask)) != 0xFFFFFFFFFFFFFFFF) {
          break;
        }
      #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        if ((mss_a | mask) != 0xFFFFFFFFFFFFFFFF) {
          break;
        }
      #else
        const int leading_ones = get_leading_ones(mss + j, 10);
        if (leading_ones < leading_zeros) {
          break;
        }
      #endif
      const uint64_t a = usr_a | mss_a;
      const uint16_t b = usr_b | mss_b;
      if (a == 0xFFFFFFFFFFFFFFFF && b == 0xFFFF) {
        return 1;
      }
    }
  }
  return 0;
}

int main() {
  const char usr_path[] = NCUBE_DIR "5_usr_2.bin";
  const char mss_path[] = NCUBE_DIR "5_mss_2.bin";
  char *usr, *mss;
  const size_t usr_len = read_from_file(usr_path, &usr);
  if (usr_len == 0) {
    return 1;
  }
  const size_t mss_len = read_from_file(mss_path, &mss);
  if (mss_len == 0) {
    return 2;
  }
  const int slices_all = combine_usr_mss(usr, usr_len, mss, mss_len);
  printf("%d\n", slices_all);
}
