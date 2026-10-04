#pragma once

#include <stdio.h>
#include <alsa/asoundlib.h>
#include "SoundChainPlatform.hpp"

class AlsaSCP : public SoundChainPlatform {
public:
	AlsaSCP() {}
	~AlsaSCP() {}

private:
	snd_pcm_t *pcm_handle;

	void Setup() override {
		snd_pcm_hw_params_t *params;
		int err;

		if ((err = snd_pcm_open(&pcm_handle, PCM_DEVICE, SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
			fprintf(stderr, "Unable to open PCM device: %s\n", snd_strerror(err));
			return;
		}

		snd_pcm_hw_params_malloc(&params);
		snd_pcm_hw_params_any(pcm_handle, params);
		
		snd_pcm_hw_params_set_access(pcm_handle, params, SND_PCM_ACCESS_RW_INTERLEAVED);
		snd_pcm_hw_params_set_format(pcm_handle, params, SND_PCM_FORMAT_FLOAT);

		unsigned int rate = GetSoundChainSettings().SampleRate;
		snd_pcm_hw_params_set_rate_near(pcm_handle, params, &rate, 0);
		// _settings.SampleRate = rate;
		snd_pcm_hw_params_set_channels(pcm_handle, params, GetSoundChainSettings().Channels);

		if ((err = snd_pcm_hw_params(pcm_handle, params)) < 0) {
			fprintf(stderr, "Unable to set hardware parameters: %s\n", snd_strerror(err));
			return;
		}
		snd_pcm_hw_params_free(params);


	}

	void Start() override {
		int err;

		if ((err = snd_pcm_prepare(pcm_handle)) < 0) {
	        fprintf(stderr, "Unable to prepare PCM device: %s\n", snd_strerror(err));
	        return;
	    }
		
		if ((err = snd_pcm_start(pcm_handle)) < 0) {
	        fprintf(stderr, "Unable to start PCM device: %s\n", snd_strerror(err));
	        return;
	    }
	}

	void End() override {
		snd_pcm_drain(pcm_handle);
		snd_pcm_close(pcm_handle);
	}
};