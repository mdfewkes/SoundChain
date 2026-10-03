#pragma once

#include <alsa/asoundlib.h>
#include "SoundChainPlatform.hpp"

class AlsaSCP : public SoundChainPlatform {
public:
	AlsaSCP() {}
	~AlsaSCP() {}

private:
	snd_pcm_t *handle;

	void Setup() override {
	}

	void Start() override {
		snd_pcm_open(&handle, "default", SND_PCM_STREAM_PLAYBACK, 0);
	}

	void End() override {
		snd_pcm_close(handle);
	}
};