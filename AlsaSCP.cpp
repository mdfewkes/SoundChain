#pragma once

#include <stdio.h>
#include <alsa/asoundlib.h>
#include <thread>
#include <atomic>
#include "SoundChainPlatform.hpp"

class AlsaSCP : public SoundChainPlatform {
public:
	AlsaSCP() {}
	~AlsaSCP() {
		if (buffer) delete[] buffer;
	}

	void AudioThread(AlsaSCP* alsaSCP) {
		if (!buffer) return;

		while (running.load()) {
			FillBuffer(buffer, periodSize);

			snd_pcm_sframes_t err = snd_pcm_writei(pcm_handle, buffer, periodSize);
			if (err == -EPIPE) {    /* under-run */
				err = snd_pcm_prepare(pcm_handle);
				if (err < 0)  {
					printf("Can't recovery from underrun, prepare failed: %s\n", snd_strerror(err));
					return;
				}
			} else if (err == -ESTRPIPE) {
				while ((err = snd_pcm_resume(pcm_handle)) == -EAGAIN) {
					usleep(250000);
				}
				if (err < 0) {
					err = snd_pcm_prepare(pcm_handle);
					if (err < 0)
						printf("Can't recovery from suspend, prepare failed: %s\n", snd_strerror(err));
				}
			}
		}
	}

private:
	const int BUFFER_PERIOD = 512;
	snd_pcm_t *pcm_handle;
	std::thread audio_thread;
	std::atomic<bool> running;
	float* buffer;
	int periodSize;

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

		snd_pcm_hw_params_set_channels(pcm_handle, params, GetSoundChainSettings().Channels);
		unsigned int rate = GetSoundChainSettings().SampleRate;
		snd_pcm_hw_params_set_rate_near(pcm_handle, params, &rate, 0);
		// _settings.SampleRate = rate;

		snd_pcm_uframes_t periodSize = (snd_pcm_uframes_t)BUFFER_PERIOD;
		snd_pcm_hw_params_set_period_time_near(pcm_handle, params, &periodSize, NULL);
		snd_pcm_uframes_t bufferSize = periodSize * 4;
		snd_pcm_hw_params_set_buffer_size_near(pcm_handle, params, &bufferSize);
		buffer = new float[bufferSize];

		if ((err = snd_pcm_hw_params(pcm_handle, params)) < 0) {
			fprintf(stderr, "Unable to set hardware parameters: %s\n", snd_strerror(err));
			return;
		}
		snd_pcm_hw_params_free(params);

		if ((err = snd_pcm_prepare(pcm_handle)) < 0) {
			fprintf(stderr, "Unable to prepare PCM device: %s\n", snd_strerror(err));
			return;
		}

	}

	void Start() override {
		if (running.load()) return;

		running.store(true);
		audio_thread = std::thread(&AlsaSCP::AudioThread, this);
	}

	void End() override {

		running.store(false);
		audio_thread.join();

		snd_pcm_drain(pcm_handle);
		snd_pcm_close(pcm_handle);
	}
};