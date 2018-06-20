#pragma once

#include <stdbool.h>

struct cr_ft_state;

struct game_api {
    /**
     * @return a fresh game state
     */
    struct cr_ft_state *(*init)(bool bInteractive);

    /**
     * Destroys a game state.
     */
    void (*finalize)(struct cr_ft_state *state);

    /**
     * Called exactly once when the game code is reloaded.
     */
    void (*reload)(struct cr_ft_state *state);

    /**
     * Called exactly once when the game code is about to be reloaded.
     */
    void (*unload)(struct cr_ft_state *state);

    /**
     * Called at a regular interval by the main program.
     * @return true if the program should continue
     */
    bool (*step)(struct cr_ft_state *state);
};

extern const struct game_api GAME_API;
// extern struct game_api GAME_API;
