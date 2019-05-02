#pragma once

// Notice that index "_with_line_ind" starts from 1:
#define WITH_LINES(_with_lines_name, _with_lines_ind, _with_lines_tokens, _with_lines_codeblock...) \
  do \
    {                                                                   \
      std::vector<std::string> _with_lines_tokens;                      \
      FILE * _with_lines_fp;                                            \
      char * _with_lines_line = NULL;                                   \
      size_t _with_lines_len = 0;                                       \
      ssize_t _with_lines_read;                                         \
      int _with_lines_ind = 0;                                          \
      _with_lines_fp = fopen(_with_lines_name, "r");                    \
      if (_with_lines_fp == NULL)                                       \
        exit(EXIT_FAILURE);                                             \
      while ((_with_lines_read = getline(&_with_lines_line, &_with_lines_len, _with_lines_fp)) != -1) { \
        _with_lines_ind++;                                              \
        _with_lines_tokens = str_split(_with_lines_line,' ');           \
        _with_lines_codeblock                                           \
          }                                                             \
      fclose(_with_lines_fp);                                           \
      if (_with_lines_line)                                             \
        free(_with_lines_line);                                         \
    }                                                                   \
  while(0)

#define SET3_IND(_SET3_IND_TARGET_ROOT_NAME, _SET3_IND_NAME_SPECIFIER1, _SET3_IND_NAME_SPECIFIER2, _SET3_IND_NAME_SPECIFIER3, _SET3_IND_TARGET_CATEGORY, _SET3_IND_SRC_IND) \
do									\
{									\
 _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER1##_SET3_IND_TARGET_CATEGORY = _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER1[_SET3_IND_SRC_IND]; \
 _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER2##_SET3_IND_TARGET_CATEGORY = _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER2[_SET3_IND_SRC_IND]; \
 _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER3##_SET3_IND_TARGET_CATEGORY = _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER3[_SET3_IND_SRC_IND]; \
 } \
while(0)
