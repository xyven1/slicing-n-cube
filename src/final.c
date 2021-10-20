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

void print_binary(const char *buf, size_t num_bytes) {
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

int combine_usr_mss(const char *usr, size_t num_usr, const char *mss,
                    size_t num_mss) {
  for (size_t i = 0; i < num_usr; i += 10) {
    for (size_t j = 0; j < num_mss; j += 10) {
      const uint64_t usr_a = *(const uint64_t *)(usr + i);
      const uint64_t mss_a = *(const uint64_t *)(mss + j);
      const uint16_t usr_b = *(const uint16_t *)(usr + i + 8);
      const uint16_t mss_b = *(const uint16_t *)(mss + j + 8);
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
  const size_t num_usr = read_from_file(usr_path, &usr);
  if (num_usr == 0) {
    return 1;
  }
  const size_t num_mss = read_from_file(mss_path, &mss);
  if (num_mss == 0) {
    return 2;
  }
  const int slices_all = combine_usr_mss(usr, num_usr, mss, num_mss);
  printf("%d\n", slices_all);
}
